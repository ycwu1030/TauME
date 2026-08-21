#ifndef TAUAMP_DELPHES_PION_PAIR_ROOT_WRITER_H_
#define TAUAMP_DELPHES_PION_PAIR_ROOT_WRITER_H_

#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "tauamp/delphes/pion_pair_event_loop.h"

namespace tauamp::delphes {

struct PionPairRootOutputMetadata {
    std::string input_file_sha256;
    std::string input_tree_name;
    PionObjectSource object_source;
    double electron_polarization;
    double positron_polarization;
    double tau_mass;
    tautau::ElectroweakParameters electroweak_parameters;
    std::string implementation_revision;
};

namespace detail {

inline void validate_manifest_value(const std::string& value, const char* name) {
    if (value.empty()) throw std::invalid_argument(std::string(name) + " must be non-empty");
    if (value.find_first_of("\r\n") != std::string::npos)
        throw std::invalid_argument(std::string(name) + " must not contain a line break");
}

inline const char* object_source_name(PionObjectSource source) {
    switch (source) {
        case PionObjectSource::generated:
            return "generated";
        case PionObjectSource::reconstructed:
            return "reconstructed";
    }
    throw std::invalid_argument("pion object source is invalid");
}

inline PionPairSelectionStatus selection_status(tautau::HadronicTauPairSolutionStatus status) {
    switch (status) {
        case tautau::HadronicTauPairSolutionStatus::no_solution:
            return PionPairSelectionStatus::kinematic_no_solution;
        case tautau::HadronicTauPairSolutionStatus::non_unique:
            return PionPairSelectionStatus::kinematic_non_unique;
        case tautau::HadronicTauPairSolutionStatus::one_solution:
            return PionPairSelectionStatus::evaluable_one_solution;
        case tautau::HadronicTauPairSolutionStatus::two_solutions:
            return PionPairSelectionStatus::evaluable_two_solutions;
    }
    throw std::invalid_argument("hadronic tau-pair solution status is invalid");
}

inline void validate_metadata(const PionPairRootOutputMetadata& metadata) {
    validate_manifest_value(metadata.input_file_sha256, "input file SHA-256");
    validate_manifest_value(metadata.input_tree_name, "input tree name");
    validate_manifest_value(metadata.implementation_revision, "implementation revision");
    if (!std::isfinite(metadata.electron_polarization) || !std::isfinite(metadata.positron_polarization) ||
        std::abs(metadata.electron_polarization) > 1.0 || std::abs(metadata.positron_polarization) > 1.0)
        throw std::invalid_argument("beam polarizations must be finite and lie in [-1, 1]");
    if (!std::isfinite(metadata.tau_mass) || metadata.tau_mass <= 0.0)
        throw std::invalid_argument("tau mass must be finite and positive");
    const tautau::ElectroweakParameters& parameters = metadata.electroweak_parameters;
    if (!std::isfinite(parameters.electron_charge) || !std::isfinite(parameters.sin_theta_w) ||
        !std::isfinite(parameters.z_mass) || !std::isfinite(parameters.z_width))
        throw std::invalid_argument("electroweak parameters must be finite");
}

inline std::string metadata_manifest(const PionPairRootOutputMetadata& metadata) {
    std::ostringstream manifest;
    manifest << std::setprecision(17);
    manifest << "schema_version=1\n";
    manifest << "input_file_sha256=" << metadata.input_file_sha256 << '\n';
    manifest << "input_tree_name=" << metadata.input_tree_name << '\n';
    manifest << "object_source=" << object_source_name(metadata.object_source) << '\n';
    manifest << "electron_polarization=" << metadata.electron_polarization << '\n';
    manifest << "positron_polarization=" << metadata.positron_polarization << '\n';
    manifest << "tau_mass=" << metadata.tau_mass << '\n';
    manifest << "electric_charge=" << metadata.electroweak_parameters.electron_charge << '\n';
    manifest << "weak_mixing_sine=" << metadata.electroweak_parameters.sin_theta_w << '\n';
    manifest << "z_mass=" << metadata.electroweak_parameters.z_mass << '\n';
    manifest << "z_width=" << metadata.electroweak_parameters.z_width << '\n';
    manifest << "momentum_frame=lab\n";
    manifest << "beam_ordering=electron,positron\n";
    manifest << "pion_ordering=pion_minus,pion_plus\n";
    manifest << "component_ordering=sm,f2_real,f2_imaginary,f3_real,f3_imaginary\n";
    manifest << "branch_weight_rule=equal_kinematic_weight\n";
    manifest << "truth_tau_input_used=false\n";
    manifest << "classification_encoding=supported_pion_pair:0,recognized_unsupported:1,ambiguous:2,unclassified:3\n";
    manifest << "tau_momentum_provenance_encoding=generator_truth_available:0,unavailable:1\n";
    manifest << "eligibility_status_encoding=unsupported_input:0,kinematic_no_solution:1,kinematic_non_unique:2,evaluable_one_solution:3,evaluable_two_solutions:4\n";
    manifest << "implementation_revision=" << metadata.implementation_revision << '\n';
    return manifest.str();
}

}  // namespace detail

class PionPairRootWriter {
public:
    PionPairRootWriter(const std::string& output_path, PionPairRootOutputMetadata metadata) : metadata_(std::move(metadata)) {
        detail::validate_metadata(metadata_);
        file_.reset(TFile::Open(output_path.c_str(), "CREATE"));
        if (!file_ || file_->IsZombie()) throw std::runtime_error("cannot create pion-pair ROOT output file");
        file_->cd();
        events_ = new TTree("Events", "Bounded pion-pair event records");
        hypotheses_ = new TTree("Hypotheses", "Bounded pion-pair kinematic hypotheses");
        bind_event_branches();
        bind_hypothesis_branches();
    }

    PionPairRootWriter(const PionPairRootWriter&) = delete;
    PionPairRootWriter& operator=(const PionPairRootWriter&) = delete;
    PionPairRootWriter(PionPairRootWriter&&) = delete;
    PionPairRootWriter& operator=(PionPairRootWriter&&) = delete;

    ~PionPairRootWriter() {
        try {
            finalize();
        } catch (...) {
        }
    }

    void write(const PionPairEventLoopRecord& record) {
        if (finalized_) throw std::logic_error("cannot write a finalized pion-pair ROOT output");
        if (has_previous_entry_ && record.entry <= previous_entry_)
            throw std::invalid_argument("pion-pair ROOT output entries must be strictly increasing");

        const Int_t hypothesis_count = validate_record(record);
        load_event(record, hypothesis_count);
        events_->Fill();
        if (record.evaluation.has_value()) {
            const auto& entries = record.evaluation->entries();
            for (std::size_t index = 0; index < entries.size(); ++index) {
                load_hypothesis(record.entry, static_cast<Int_t>(index), entries[index]);
                hypotheses_->Fill();
            }
        }
        previous_entry_ = record.entry;
        has_previous_entry_ = true;
    }

    void finalize() {
        if (finalized_) return;
        file_->cd();
        TObjString manifest(detail::metadata_manifest(metadata_).c_str());
        manifest.Write("TauAmpPionPairMetadata");
        events_->Write();
        hypotheses_->Write();
        file_->Write();
        file_->Close();
        finalized_ = true;
    }

private:
    void bind_event_branches() {
        events_->Branch("source_entry", &event_source_entry_);
        events_->Branch("classification", &event_classification_);
        events_->Branch("tau_momentum_provenance", &event_tau_momentum_provenance_);
        events_->Branch("eligibility_status", &event_eligibility_status_);
        events_->Branch("hypothesis_count", &event_hypothesis_count_);
        events_->Branch("has_lab_observation", &event_has_lab_observation_);
        events_->Branch("electron_px", &electron_px_);
        events_->Branch("electron_py", &electron_py_);
        events_->Branch("electron_pz", &electron_pz_);
        events_->Branch("electron_e", &electron_e_);
        events_->Branch("positron_px", &positron_px_);
        events_->Branch("positron_py", &positron_py_);
        events_->Branch("positron_pz", &positron_pz_);
        events_->Branch("positron_e", &positron_e_);
        events_->Branch("pion_minus_px", &pion_minus_px_);
        events_->Branch("pion_minus_py", &pion_minus_py_);
        events_->Branch("pion_minus_pz", &pion_minus_pz_);
        events_->Branch("pion_minus_e", &pion_minus_e_);
        events_->Branch("pion_plus_px", &pion_plus_px_);
        events_->Branch("pion_plus_py", &pion_plus_py_);
        events_->Branch("pion_plus_pz", &pion_plus_pz_);
        events_->Branch("pion_plus_e", &pion_plus_e_);
    }

    void bind_hypothesis_branches() {
        hypotheses_->Branch("source_entry", &hypothesis_source_entry_);
        hypotheses_->Branch("hypothesis_index", &hypothesis_index_);
        hypotheses_->Branch("kinematic_weight", &kinematic_weight_);
        hypotheses_->Branch("component_sm", &component_sm_);
        hypotheses_->Branch("component_f2_real", &component_f2_real_);
        hypotheses_->Branch("component_f2_imaginary", &component_f2_imaginary_);
        hypotheses_->Branch("component_f3_real", &component_f3_real_);
        hypotheses_->Branch("component_f3_imaginary", &component_f3_imaginary_);
    }

    static Int_t validate_record(const PionPairEventLoopRecord& record) {
        if (record.classification != TauPairDecayClassification::supported_pion_pair) {
            if (record.selection_status != PionPairSelectionStatus::unsupported_input || record.evaluation.has_value() ||
                record.observation.has_value())
                throw std::invalid_argument("unsupported pion-pair record is internally inconsistent");
            return 0;
        }
        if (!record.observation.has_value() || !record.evaluation.has_value())
            throw std::invalid_argument("supported pion-pair record lacks a lab observation or evaluation");
        const auto& evaluation = *record.evaluation;
        if (record.selection_status != detail::selection_status(evaluation.status()))
            throw std::invalid_argument("pion-pair record eligibility status disagrees with kinematics");
        const std::size_t count = evaluation.entries().size();
        const std::size_t expected = evaluation.status() == tautau::HadronicTauPairSolutionStatus::one_solution   ? 1U
                                     : evaluation.status() == tautau::HadronicTauPairSolutionStatus::two_solutions ? 2U
                                                                                                                       : 0U;
        if (count != expected) throw std::invalid_argument("pion-pair record has incompatible hypothesis count");
        return static_cast<Int_t>(count);
    }

    static void copy_momentum(const tautau::FourMomentum& source, Double_t& px, Double_t& py, Double_t& pz,
                              Double_t& energy) {
        px = source.px();
        py = source.py();
        pz = source.pz();
        energy = source.energy();
    }

    void load_event(const PionPairEventLoopRecord& record, Int_t hypothesis_count) {
        event_source_entry_ = record.entry;
        event_classification_ = static_cast<Int_t>(record.classification);
        event_tau_momentum_provenance_ = static_cast<Int_t>(record.tau_momentum_provenance);
        event_eligibility_status_ = static_cast<Int_t>(record.selection_status);
        event_hypothesis_count_ = hypothesis_count;
        event_has_lab_observation_ = record.observation.has_value();
        if (!record.observation.has_value()) {
            const Double_t nan = std::numeric_limits<Double_t>::quiet_NaN();
            electron_px_ = electron_py_ = electron_pz_ = electron_e_ = nan;
            positron_px_ = positron_py_ = positron_pz_ = positron_e_ = nan;
            pion_minus_px_ = pion_minus_py_ = pion_minus_pz_ = pion_minus_e_ = nan;
            pion_plus_px_ = pion_plus_py_ = pion_plus_pz_ = pion_plus_e_ = nan;
            return;
        }
        const PionPairLabObservation& observation = *record.observation;
        copy_momentum(observation.beams.electron_lab(), electron_px_, electron_py_, electron_pz_, electron_e_);
        copy_momentum(observation.beams.positron_lab(), positron_px_, positron_py_, positron_pz_, positron_e_);
        copy_momentum(observation.pion_minus_lab, pion_minus_px_, pion_minus_py_, pion_minus_pz_, pion_minus_e_);
        copy_momentum(observation.pion_plus_lab, pion_plus_px_, pion_plus_py_, pion_plus_pz_, pion_plus_e_);
    }

    void load_hypothesis(Long64_t source_entry, Int_t index, const tautau::PionPairHypothesisComponents& entry) {
        hypothesis_source_entry_ = source_entry;
        hypothesis_index_ = index;
        kinematic_weight_ = entry.weight;
        component_sm_ = entry.components.sm;
        component_f2_real_ = entry.components.f2_real;
        component_f2_imaginary_ = entry.components.f2_imaginary;
        component_f3_real_ = entry.components.f3_real;
        component_f3_imaginary_ = entry.components.f3_imaginary;
    }

    PionPairRootOutputMetadata metadata_;
    std::unique_ptr<TFile> file_;
    TTree* events_{};
    TTree* hypotheses_{};
    bool finalized_{};
    bool has_previous_entry_{};
    Long64_t previous_entry_{};

    Long64_t event_source_entry_{};
    Int_t event_classification_{};
    Int_t event_tau_momentum_provenance_{};
    Int_t event_eligibility_status_{};
    Int_t event_hypothesis_count_{};
    Bool_t event_has_lab_observation_{};
    Double_t electron_px_{};
    Double_t electron_py_{};
    Double_t electron_pz_{};
    Double_t electron_e_{};
    Double_t positron_px_{};
    Double_t positron_py_{};
    Double_t positron_pz_{};
    Double_t positron_e_{};
    Double_t pion_minus_px_{};
    Double_t pion_minus_py_{};
    Double_t pion_minus_pz_{};
    Double_t pion_minus_e_{};
    Double_t pion_plus_px_{};
    Double_t pion_plus_py_{};
    Double_t pion_plus_pz_{};
    Double_t pion_plus_e_{};

    Long64_t hypothesis_source_entry_{};
    Int_t hypothesis_index_{};
    Double_t kinematic_weight_{};
    Double_t component_sm_{};
    Double_t component_f2_real_{};
    Double_t component_f2_imaginary_{};
    Double_t component_f3_real_{};
    Double_t component_f3_imaginary_{};
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_PION_PAIR_ROOT_WRITER_H_
