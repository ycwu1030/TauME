#ifndef TAUAMP_DELPHES_COMMON_CHANNEL_DELPHES_READER_H_
#define TAUAMP_DELPHES_COMMON_CHANNEL_DELPHES_READER_H_

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include "classes/DelphesClasses.h"
#include "tauamp/tautau/ordered_channel_adapter.h"

namespace tauamp::delphes {

struct CommonChannelDelphesReaderConfig {
    tautau::FourMomentum beam_electron;
    tautau::FourMomentum beam_positron;
};

struct CommonChannelPionCollections {
    std::vector<tautau::ReconstructedPion> charged;
    std::vector<tautau::ReconstructedPion> neutral;
};

class CommonChannelDelphesReader {
public:
    CommonChannelDelphesReader(const std::string& file_path, CommonChannelDelphesReaderConfig config)
        : config_(config) {
        file_.reset(TFile::Open(file_path.c_str(), "READ"));
        if (!file_ || file_->IsZombie()) throw std::runtime_error("cannot open common-channel Delphes ROOT file");
        tree_ = dynamic_cast<TTree*>(file_->Get("Delphes"));
        if (tree_ == nullptr) throw std::runtime_error("common-channel Delphes ROOT file has no Delphes tree");
        bind_branch("ChargedPion", charged_pions_);
        bind_branch("NeutralPion", neutral_pions_);
    }

    Long64_t entry_count() const { return tree_->GetEntries(); }

    CommonChannelPionCollections read(Long64_t entry) const {
        if (entry < 0 || entry >= entry_count()) throw std::out_of_range("Delphes entry is outside tree range");
        if (tree_->GetEntry(entry) < 0) throw std::runtime_error("cannot read Delphes tree entry");
        return {read_collection(*charged_pions_, "ChargedPion"), read_collection(*neutral_pions_, "NeutralPion")};
    }

private:
    void bind_branch(const char* name, TClonesArray*& target) {
        if (tree_->GetBranch(name) == nullptr) throw std::runtime_error(std::string("Delphes tree is missing branch ") + name);
        if (tree_->SetBranchAddress(name, &target) < 0)
            throw std::runtime_error(std::string("cannot bind Delphes branch ") + name);
    }

    static tautau::ReconstructedPion convert(const GenParticle& particle, std::size_t index) {
        if (particle.PID != 211 && particle.PID != -211 && particle.PID != 111)
            throw std::invalid_argument("common-channel pion collection contains a non-pion PID");
        const double pt = std::hypot(particle.Px, particle.Py);
        const double momentum = std::sqrt(particle.Px * particle.Px + particle.Py * particle.Py + particle.Pz * particle.Pz);
        const double eta = momentum == std::abs(particle.Pz) ? (particle.Pz >= 0.0 ? 1.0e9 : -1.0e9)
                                                              : 0.5 * std::log((momentum + particle.Pz) / (momentum - particle.Pz));
        const double phi = std::atan2(particle.Py, particle.Px);
        return {index, particle.PID, pt, eta, phi, {particle.Px, particle.Py, particle.Pz, particle.E}};
    }

    static std::vector<tautau::ReconstructedPion> read_collection(const TClonesArray& collection, const char* name) {
        std::vector<tautau::ReconstructedPion> result;
        result.reserve(static_cast<std::size_t>(collection.GetEntriesFast()));
        for (int index = 0; index < collection.GetEntriesFast(); ++index) {
            const auto* particle = dynamic_cast<const GenParticle*>(collection.At(index));
            if (particle == nullptr) throw std::runtime_error(std::string("Delphes ") + name + " contains a non-GenParticle object");
            result.push_back(convert(*particle, static_cast<std::size_t>(index)));
        }
        return result;
    }

    CommonChannelDelphesReaderConfig config_;
    std::unique_ptr<TFile> file_;
    TTree* tree_{nullptr};
    TClonesArray* charged_pions_{nullptr};
    TClonesArray* neutral_pions_{nullptr};
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_COMMON_CHANNEL_DELPHES_READER_H_
