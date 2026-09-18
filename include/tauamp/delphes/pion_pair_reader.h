#ifndef TAUAMP_DELPHES_PION_PAIR_READER_H_
#define TAUAMP_DELPHES_PION_PAIR_READER_H_

#include <cstddef>
#include <cmath>
#include <initializer_list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include "classes/DelphesClasses.h"
#include "tauamp/tautau/kinematics.h"

namespace tauamp::delphes {

enum class PionObjectSource { generated, reconstructed, truth, truth_neutrino, truth_tau_reconstructed_pions, truth_neutrino_reconstructed_pions };

enum class TauPairDecayClassification { supported_pion_pair, recognized_unsupported, ambiguous, unclassified };

enum class TauMomentumProvenance { generator_truth_available, unavailable };

struct PionCandidate {
    int pid;
    int charge;
};

inline bool is_physical_charged_pion(const PionCandidate& candidate) {
    return (candidate.pid == -211 && candidate.charge == -1) || (candidate.pid == 211 && candidate.charge == 1);
}

inline TauPairDecayClassification classify_pion_pair(std::initializer_list<PionCandidate> candidates,
                                                     std::size_t neutral_pion_count) {
    std::size_t pion_minus_count = 0U;
    std::size_t pion_plus_count = 0U;
    for (const PionCandidate& candidate : candidates) {
        if (!is_physical_charged_pion(candidate)) return TauPairDecayClassification::unclassified;
        if (candidate.pid == -211) {
            ++pion_minus_count;
        } else {
            ++pion_plus_count;
        }
    }
    if (pion_minus_count == 1U && pion_plus_count == 1U && candidates.size() == 2U)
        return neutral_pion_count == 0U ? TauPairDecayClassification::supported_pion_pair
                                        : TauPairDecayClassification::recognized_unsupported;
    if (pion_minus_count > 0U && pion_plus_count > 0U) return TauPairDecayClassification::ambiguous;
    return TauPairDecayClassification::unclassified;
}

struct PionPairReaderConfig {
    double electron_polarization;
    double positron_polarization;
    double tau_mass;
    PionObjectSource source;
};

struct PionPairLabObservation {
    tautau::BeamState beams;
    tautau::FourMomentum pion_minus_lab;
    tautau::FourMomentum pion_plus_lab;
    double tau_mass;
    std::optional<tautau::FourMomentum> truth_tau_minus_lab;
    std::optional<tautau::FourMomentum> truth_tau_plus_lab;
};

struct PionPairDelphesEvent {
    TauPairDecayClassification classification;
    TauMomentumProvenance tau_momentum_provenance;
    std::optional<PionPairLabObservation> observation;
};

class PionPairDelphesReader {
public:
    PionPairDelphesReader(const std::string& file_path, PionPairReaderConfig config) : config_(config) {
        if (config_.tau_mass <= 0.0) throw std::invalid_argument("tau mass must be positive");
        if ((config_.source == PionObjectSource::truth || config_.source == PionObjectSource::truth_neutrino || config_.source == PionObjectSource::truth_tau_reconstructed_pions || config_.source == PionObjectSource::truth_neutrino_reconstructed_pions) && std::abs(config_.tau_mass - 1.777) > 1e-12)
            throw std::invalid_argument("truth input route requires MG5/Delphes tau mass 1.777 GeV");
        file_.reset(TFile::Open(file_path.c_str(), "READ"));
        if (!file_ || file_->IsZombie()) throw std::runtime_error("cannot open Delphes ROOT file");
        tree_ = dynamic_cast<TTree*>(file_->Get("Delphes"));
        if (tree_ == nullptr) throw std::runtime_error("Delphes ROOT file has no Delphes tree");
        bind_branch("BeamParticle", beam_particles_);
        bind_branch("GenTau", generator_taus_);
        bind_branch("GenNeutrino", generator_neutrinos_);
        bind_branch("GenChargedPion", generator_charged_pions_);
        bind_branch("GenNeutralPion", generator_neutral_pions_);
        bind_branch("ChargedPion", reconstructed_charged_pions_);
        bind_branch("NeutralPion", reconstructed_neutral_pions_);
    }

    PionPairDelphesReader(const PionPairDelphesReader&) = delete;
    PionPairDelphesReader& operator=(const PionPairDelphesReader&) = delete;
    PionPairDelphesReader(PionPairDelphesReader&&) = delete;
    PionPairDelphesReader& operator=(PionPairDelphesReader&&) = delete;

    Long64_t entry_count() const { return tree_->GetEntries(); }

    PionPairDelphesEvent read(Long64_t entry) const {
        if (entry < 0 || entry >= tree_->GetEntries()) throw std::out_of_range("Delphes entry is outside tree range");
        if (tree_->GetEntry(entry) < 0) throw std::runtime_error("cannot read Delphes tree entry");

        const tautau::BeamState beams = read_beams();
        const bool use_generated_pions = config_.source == PionObjectSource::generated || config_.source == PionObjectSource::truth || config_.source == PionObjectSource::truth_neutrino;
        const TClonesArray* charged_pions = use_generated_pions ? generator_charged_pions_ : reconstructed_charged_pions_;
        const TClonesArray* neutral_pions = use_generated_pions ? generator_neutral_pions_ : reconstructed_neutral_pions_;
        const std::vector<const GenParticle*> charged = objects(*charged_pions, "charged pion");
        const TauPairDecayClassification classification =
            classify(charged, static_cast<std::size_t>(neutral_pions->GetEntriesFast()));
        const TauMomentumProvenance provenance =
            config_.source == PionObjectSource::reconstructed ? TauMomentumProvenance::unavailable : generator_truth_provenance();

        if (classification != TauPairDecayClassification::supported_pion_pair)
            return {classification, provenance, std::nullopt};
        if (config_.source == PionObjectSource::truth || config_.source == PionObjectSource::truth_tau_reconstructed_pions)
            return {classification, provenance, make_truth_observation(beams, charged)};
        if (config_.source == PionObjectSource::truth_neutrino || config_.source == PionObjectSource::truth_neutrino_reconstructed_pions)
            return {classification, provenance, make_neutrino_observation(beams, charged)};
        return {classification, provenance, make_observation(beams, charged)};
    }

private:
    void bind_branch(const char* name, TClonesArray*& target) {
        if (tree_->GetBranch(name) == nullptr) throw std::runtime_error(std::string("Delphes tree is missing branch ") + name);
        if (tree_->SetBranchAddress(name, &target) < 0)
            throw std::runtime_error(std::string("cannot bind Delphes branch ") + name);
    }

    static const GenParticle& particle_at(const TClonesArray& collection, int index, const char* collection_name) {
        const auto* particle = dynamic_cast<const GenParticle*>(collection.At(index));
        if (particle == nullptr) throw std::runtime_error(std::string("Delphes ") + collection_name + " contains a non-GenParticle object");
        return *particle;
    }

    static std::vector<const GenParticle*> objects(const TClonesArray& collection, const char* collection_name) {
        std::vector<const GenParticle*> result;
        const int count = collection.GetEntriesFast();
        result.reserve(static_cast<std::size_t>(count));
        for (int index = 0; index < count; ++index) result.push_back(&particle_at(collection, index, collection_name));
        return result;
    }

    static tautau::FourMomentum momentum(const GenParticle& particle) {
        return {particle.Px, particle.Py, particle.Pz, particle.E};
    }

    tautau::BeamState read_beams() const {
        const std::vector<const GenParticle*> candidates = objects(*beam_particles_, "BeamParticle");
        const GenParticle* electron = nullptr;
        const GenParticle* positron = nullptr;
        for (const GenParticle* candidate : candidates) {
            if (candidate->PID == 11) {
                if (electron != nullptr) throw std::runtime_error("BeamParticle has multiple electrons");
                electron = candidate;
            } else if (candidate->PID == -11) {
                if (positron != nullptr) throw std::runtime_error("BeamParticle has multiple positrons");
                positron = candidate;
            }
        }
        if (electron == nullptr || positron == nullptr) throw std::runtime_error("BeamParticle lacks an electron-positron pair");
        return {momentum(*electron), momentum(*positron), config_.electron_polarization, config_.positron_polarization};
    }

    static TauPairDecayClassification classify(const std::vector<const GenParticle*>& charged_pions,
                                               std::size_t neutral_pion_count) {
        std::vector<PionCandidate> candidates;
        candidates.reserve(charged_pions.size());
        for (const GenParticle* pion : charged_pions) candidates.push_back({pion->PID, pion->Charge});
        std::size_t pion_minus_count = 0U;
        std::size_t pion_plus_count = 0U;
        for (const PionCandidate& candidate : candidates) {
            if (!is_physical_charged_pion(candidate)) return TauPairDecayClassification::unclassified;
            if (candidate.pid == -211) {
                ++pion_minus_count;
            } else {
                ++pion_plus_count;
            }
        }
        if (pion_minus_count == 1U && pion_plus_count == 1U && candidates.size() == 2U)
            return neutral_pion_count == 0U ? TauPairDecayClassification::supported_pion_pair
                                            : TauPairDecayClassification::recognized_unsupported;
        if (pion_minus_count > 0U && pion_plus_count > 0U) return TauPairDecayClassification::ambiguous;
        return TauPairDecayClassification::unclassified;
    }

    TauMomentumProvenance generator_truth_provenance() const {
        const std::vector<const GenParticle*> candidates = objects(*generator_taus_, "GenTau");
        const GenParticle* tau_minus = nullptr;
        const GenParticle* tau_plus = nullptr;
        for (const GenParticle* candidate : candidates) {
            if (candidate->PID == 15) {
                if (tau_minus != nullptr) throw std::runtime_error("GenTau has multiple tau-minus objects");
                tau_minus = candidate;
            } else if (candidate->PID == -15) {
                if (tau_plus != nullptr) throw std::runtime_error("GenTau has multiple tau-plus objects");
                tau_plus = candidate;
            }
        }
        if (tau_minus == nullptr || tau_plus == nullptr) throw std::runtime_error("GenTau lacks a tau pair");
        return TauMomentumProvenance::generator_truth_available;
    }

    PionPairLabObservation make_observation(const tautau::BeamState& beams,
                                            const std::vector<const GenParticle*>& charged_pions) const {
        const GenParticle* pion_minus = nullptr;
        const GenParticle* pion_plus = nullptr;
        for (const GenParticle* pion : charged_pions) {
            if (pion->PID == -211) {
                pion_minus = pion;
            } else if (pion->PID == 211) {
                pion_plus = pion;
            }
        }
        if (pion_minus == nullptr || pion_plus == nullptr)
            throw std::logic_error("supported pion-pair classification lacks charge-labelled pion objects");
        return {beams, momentum(*pion_minus), momentum(*pion_plus), config_.tau_mass, std::nullopt, std::nullopt};
    }

    PionPairLabObservation make_truth_observation(const tautau::BeamState& beams,
                                                  const std::vector<const GenParticle*>& charged_pions) const {
        const auto pion_observation = make_observation(beams, charged_pions);
        const std::vector<const GenParticle*> taus = objects(*generator_taus_, "GenTau");
        const GenParticle* tau_minus = nullptr;
        const GenParticle* tau_plus = nullptr;
        for (const GenParticle* tau : taus) {
            if (tau->PID == 15) tau_minus = tau;
            if (tau->PID == -15) tau_plus = tau;
        }
        if (tau_minus == nullptr || tau_plus == nullptr) throw std::runtime_error("GenTau lacks a tau pair");
        const tautau::FourMomentum tau_minus_momentum = momentum(*tau_minus);
        const tautau::FourMomentum tau_plus_momentum = momentum(*tau_plus);
        const double minus_mass_squared = tau_minus_momentum.mass_squared();
        const double plus_mass_squared = tau_plus_momentum.mass_squared();
        if (minus_mass_squared <= 0.0 || plus_mass_squared <= 0.0)
            throw std::runtime_error("GenTau contains a non-time-like tau four-momentum");
        const double event_local_mass = 0.5 * (std::sqrt(minus_mass_squared) + std::sqrt(plus_mass_squared));
        return {pion_observation.beams, pion_observation.pion_minus_lab, pion_observation.pion_plus_lab,
                event_local_mass, tau_minus_momentum, tau_plus_momentum};
    }

    PionPairLabObservation make_neutrino_observation(const tautau::BeamState& beams,
                                                     const std::vector<const GenParticle*>& charged_pions) const {
        auto observation = make_observation(beams, charged_pions);
        const auto neutrinos = objects(*generator_neutrinos_, "GenNeutrino");
        const GenParticle* nu_tau = nullptr;
        const GenParticle* anti_nu_tau = nullptr;
        for (const GenParticle* neutrino : neutrinos) {
            if (neutrino->PID == 16) nu_tau = neutrino;
            if (neutrino->PID == -16) anti_nu_tau = neutrino;
        }
        if (nu_tau == nullptr || anti_nu_tau == nullptr) throw std::runtime_error("GenNeutrino lacks a tau-neutrino pair");
        const tautau::FourMomentum pion_minus = observation.pion_minus_lab;
        const tautau::FourMomentum pion_plus = observation.pion_plus_lab;
        const tautau::FourMomentum tau_minus = pion_minus + momentum(*nu_tau);
        const tautau::FourMomentum tau_plus = pion_plus + momentum(*anti_nu_tau);
        const double minus_mass_squared = tau_minus.mass_squared();
        const double plus_mass_squared = tau_plus.mass_squared();
        if (minus_mass_squared <= 0.0 || plus_mass_squared <= 0.0) throw std::runtime_error("GenNeutrino reconstruction produced non-time-like tau");
        return {observation.beams, pion_minus, pion_plus, 0.5 * (std::sqrt(minus_mass_squared) + std::sqrt(plus_mass_squared)), tau_minus, tau_plus};
    }

    PionPairReaderConfig config_;
    std::unique_ptr<TFile> file_;
    TTree* tree_{};
    TClonesArray* beam_particles_{};
    TClonesArray* generator_taus_{};
    TClonesArray* generator_neutrinos_{};
    TClonesArray* generator_charged_pions_{};
    TClonesArray* generator_neutral_pions_{};
    TClonesArray* reconstructed_charged_pions_{};
    TClonesArray* reconstructed_neutral_pions_{};
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_PION_PAIR_READER_H_
