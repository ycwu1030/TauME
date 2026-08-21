#ifndef TAUAMP_DELPHES_PION_PAIR_READER_H_
#define TAUAMP_DELPHES_PION_PAIR_READER_H_

#include <cstddef>
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

enum class PionObjectSource { generated, reconstructed };

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
        file_.reset(TFile::Open(file_path.c_str(), "READ"));
        if (!file_ || file_->IsZombie()) throw std::runtime_error("cannot open Delphes ROOT file");
        tree_ = dynamic_cast<TTree*>(file_->Get("Delphes"));
        if (tree_ == nullptr) throw std::runtime_error("Delphes ROOT file has no Delphes tree");
        bind_branch("BeamParticle", beam_particles_);
        bind_branch("GenTau", generator_taus_);
        bind_branch("GenChargedPion", generator_charged_pions_);
        bind_branch("GenNeutralPion", generator_neutral_pions_);
        bind_branch("ChargedPion", reconstructed_charged_pions_);
        bind_branch("NeutralPion", reconstructed_neutral_pions_);
    }

    PionPairDelphesReader(const PionPairDelphesReader&) = delete;
    PionPairDelphesReader& operator=(const PionPairDelphesReader&) = delete;
    PionPairDelphesReader(PionPairDelphesReader&&) = delete;
    PionPairDelphesReader& operator=(PionPairDelphesReader&&) = delete;

    PionPairDelphesEvent read(Long64_t entry) const {
        if (entry < 0 || entry >= tree_->GetEntries()) throw std::out_of_range("Delphes entry is outside tree range");
        if (tree_->GetEntry(entry) < 0) throw std::runtime_error("cannot read Delphes tree entry");

        const tautau::BeamState beams = read_beams();
        const TClonesArray* charged_pions =
            config_.source == PionObjectSource::generated ? generator_charged_pions_ : reconstructed_charged_pions_;
        const TClonesArray* neutral_pions =
            config_.source == PionObjectSource::generated ? generator_neutral_pions_ : reconstructed_neutral_pions_;
        const std::vector<const GenParticle*> charged = objects(*charged_pions, "charged pion");
        const TauPairDecayClassification classification =
            classify(charged, static_cast<std::size_t>(neutral_pions->GetEntriesFast()));
        const TauMomentumProvenance provenance =
            config_.source == PionObjectSource::generated ? generator_truth_provenance() : TauMomentumProvenance::unavailable;

        if (classification != TauPairDecayClassification::supported_pion_pair)
            return {classification, provenance, std::nullopt};
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
        return {beams, momentum(*pion_minus), momentum(*pion_plus), config_.tau_mass};
    }

    PionPairReaderConfig config_;
    std::unique_ptr<TFile> file_;
    TTree* tree_{};
    TClonesArray* beam_particles_{};
    TClonesArray* generator_taus_{};
    TClonesArray* generator_charged_pions_{};
    TClonesArray* generator_neutral_pions_{};
    TClonesArray* reconstructed_charged_pions_{};
    TClonesArray* reconstructed_neutral_pions_{};
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_PION_PAIR_READER_H_
