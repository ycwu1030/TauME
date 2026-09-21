#ifndef TAUAMP_DELPHES_COMMON_CHANNEL_ROOT_WRITER_H_
#define TAUAMP_DELPHES_COMMON_CHANNEL_ROOT_WRITER_H_

#include <array>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "tauamp/tautau/common_channel_export.h"

namespace tauamp::delphes {

class CommonChannelRootWriter {
public:
    CommonChannelRootWriter(const std::string& output_path, tautau::CommonChannelMetadata metadata)
        : metadata_(std::move(metadata)) {
        tautau::validate_common_channel_metadata(metadata_);
        file_.reset(TFile::Open(output_path.c_str(), "CREATE"));
        if (!file_ || file_->IsZombie()) throw std::runtime_error("cannot create common-channel ROOT output file");
        file_->cd();
        events_ = new TTree("Events", "G6 common ordered-channel event records");
        hypotheses_ = new TTree("Hypotheses", "G6 common ordered-channel schema-2 hypotheses");
        bind_branches();
    }

    CommonChannelRootWriter(const CommonChannelRootWriter&) = delete;
    CommonChannelRootWriter& operator=(const CommonChannelRootWriter&) = delete;

    ~CommonChannelRootWriter() {
        try {
            finalize();
        } catch (...) {
        }
    }

    void write(const tautau::CommonEventRecord& record) {
        if (finalized_) throw std::logic_error("cannot write a finalized common-channel ROOT output");
        if (has_previous_event_ && record.source_event_id <= previous_event_id_)
            throw std::invalid_argument("common-channel ROOT source event IDs must be strictly increasing");
        tautau::validate_common_event_record(metadata_, record);
        load_event(record);
        events_->Fill();
        ++event_count_;
        for (const auto& hypothesis : record.hypotheses) {
            if (!hypothesis.valid) continue;
            load_hypothesis(record.source_event_id, hypothesis);
            ++hypothesis_count_;
            hypotheses_->Fill();
        }
        previous_event_id_ = record.source_event_id;
        has_previous_event_ = true;
    }

    void finalize() {
        if (finalized_) return;
        file_->cd();
        metadata_.event_count = event_count_;
        metadata_.hypothesis_count = hypothesis_count_;
        TObjString manifest(metadata_manifest(metadata_).c_str());
        manifest.Write("TauAmpCommonChannelMetadata");
        events_->Write();
        hypotheses_->Write();
        file_->Close();
        finalized_ = true;
    }

private:
    static std::string metadata_manifest(const tautau::CommonChannelMetadata& metadata) {
        std::ostringstream output;
        output << "schema_version=" << metadata.schema_version << '\n';
        output << "model_interpretation=" << metadata.model_interpretation << '\n';
        output << "quadratic_cross_term_convention=" << metadata.quadratic_cross_term_convention << '\n';
        output << "production_bosons=" << metadata.production_bosons << '\n';
        output << "sqrt_s_GeV=" << metadata.sqrt_s_GeV << '\n';
        output << "momentum_frame=" << metadata.momentum_frame << '\n';
        output << "object_source=" << metadata.object_source << '\n';
        output << "branch_weight_rule=" << metadata.branch_weight_rule << '\n';
        output << "truth_tau_input_used=" << (metadata.truth_tau_input ? "true" : "false") << '\n';
        output << "truth_neutrino_input_used=" << (metadata.truth_neutrino_input ? "true" : "false") << '\n';
        output << "ordered_channel=" << tautau::ordered_channel_name(metadata.ordered_channel) << '\n';
        output << "tau_plus_decay=" << metadata.tau_plus_decay << '\n';
        output << "tau_minus_decay=" << metadata.tau_minus_decay << '\n';
        output << "decay_model_revision=" << metadata.decay_model_revision << '\n';
        output << "current_revision=" << metadata.current_revision << '\n';
        output << "lineshape_revision=" << metadata.lineshape_revision << '\n';
        output << "phase_space_measure_revision=" << metadata.phase_space_measure_revision << '\n';
        output << "current_normalization_revision=" << metadata.current_normalization_revision << '\n';
        output << "identical_slot_rule=" << metadata.identical_slot_rule << '\n';
        output << "charge_convention_revision=" << metadata.charge_convention_revision << '\n';
        output << "analyser_convention_revision=" << metadata.analyser_convention_revision << '\n';
        output << "pi0_object_source=" << metadata.pi0_object_source << '\n';
        output << "selection_revision=" << metadata.selection_revision << '\n';
        output << "tauoo_revision=" << metadata.tauoo_revision << '\n';
        output << "exporter_revision=" << metadata.exporter_revision << '\n';
        output << "cross_feed_policy=" << metadata.cross_feed_policy << '\n';
        output << "event_count=" << metadata.event_count << '\n';
        output << "hypothesis_count=" << metadata.hypothesis_count << '\n';
        output << "component_ordering=sm,f2_real,f2_imaginary,f3_real,f3_imaginary,f2_real_f2_real,f2_real_f2_imaginary,f2_real_f3_real,f2_real_f3_imaginary,f2_imaginary_f2_imaginary,f2_imaginary_f3_real,f2_imaginary_f3_imaginary,f3_real_f3_real,f3_real_f3_imaginary,f3_imaginary_f3_imaginary\n";
        for (const auto& [name, hash] : metadata.source_hashes) output << "source_hash." << name << '=' << hash << '\n';
        if (!metadata.failure_record.empty()) output << "failure_record=" << metadata.failure_record << '\n';
        return output.str();
    }

    void bind_branches() {
        events_->Branch("source_event_id", &event_id_);
        events_->Branch("ordered_channel", &event_channel_);
        events_->Branch("category_code", &event_category_);
        events_->Branch("truth_ordered_channel", &event_truth_channel_);
        events_->Branch("candidate_hypothesis_count", &event_hypothesis_count_);
        hypotheses_->Branch("source_event_id", &hypothesis_event_id_);
        hypotheses_->Branch("ordered_channel", &hypothesis_channel_);
        hypotheses_->Branch("hypothesis_index", &hypothesis_index_);
        hypotheses_->Branch("kinematic_weight", &hypothesis_weight_);
        hypotheses_->Branch("valid_branch_mask", &hypothesis_valid_);
        hypotheses_->Branch("components", components_.data(), "components[15]/D");
    }

    void load_event(const tautau::CommonEventRecord& record) {
        event_id_ = record.source_event_id;
        event_channel_ = static_cast<Int_t>(record.ordered_channel);
        event_category_ = record.category_code;
        event_truth_channel_ = record.truth_ordered_channel;
        event_hypothesis_count_ = static_cast<Int_t>(record.candidate_hypothesis_count);
    }

    void load_hypothesis(long long event_id, const tautau::CommonHypothesisRecord& hypothesis) {
        hypothesis_event_id_ = event_id;
        hypothesis_channel_ = static_cast<Int_t>(metadata_.ordered_channel);
        hypothesis_index_ = static_cast<Int_t>(hypothesis.hypothesis_index);
        hypothesis_weight_ = hypothesis.weight;
        hypothesis_valid_ = hypothesis.valid;
        components_ = hypothesis.components;
    }

    tautau::CommonChannelMetadata metadata_;
    std::unique_ptr<TFile> file_;
    TTree* events_{nullptr};
    TTree* hypotheses_{nullptr};
    bool finalized_{false};
    bool has_previous_event_{false};
    long long previous_event_id_{0};
    Long64_t event_id_{0};
    Int_t event_channel_{0};
    std::string event_category_;
    std::string event_truth_channel_;
    Int_t event_hypothesis_count_{0};
    Long64_t hypothesis_event_id_{0};
    Int_t hypothesis_channel_{0};
    Int_t hypothesis_index_{0};
    Double_t hypothesis_weight_{0.0};
    Bool_t hypothesis_valid_{false};
    std::array<Double_t, tautau::g6_component_count> components_{};
    long long event_count_{0};
    long long hypothesis_count_{0};
};

}  // namespace tauamp::delphes

#endif  // TAUAMP_DELPHES_COMMON_CHANNEL_ROOT_WRITER_H_
