#ifndef TAUAMP_TESTS_DELPHES_PION_PAIR_DELPHES_FIXTURE_H_
#define TAUAMP_TESTS_DELPHES_PION_PAIR_DELPHES_FIXTURE_H_

#include <memory>
#include <stdexcept>
#include <string>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include "classes/DelphesClasses.h"

namespace tauamp::delphes::test {

inline void append_fixture_particle(TClonesArray& collection, int index, int pid, int charge, double px, double py,
                                    double pz, double energy) {
    auto* particle = new (collection[index]) GenParticle();
    particle->PID = pid;
    particle->Charge = charge;
    particle->Px = static_cast<Float_t>(px);
    particle->Py = static_cast<Float_t>(py);
    particle->Pz = static_cast<Float_t>(pz);
    particle->E = static_cast<Float_t>(energy);
}

inline void write_pion_pair_delphes_fixture(const std::string& path) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
    if (!file || file->IsZombie()) throw std::runtime_error("cannot create portable Delphes fixture");

    TTree tree("Delphes", "portable pion-pair fixture");
    TClonesArray* beam_particles = new TClonesArray("GenParticle");
    TClonesArray* generator_taus = new TClonesArray("GenParticle");
    TClonesArray* generator_charged_pions = new TClonesArray("GenParticle");
    TClonesArray* generator_neutral_pions = new TClonesArray("GenParticle");
    TClonesArray* reconstructed_charged_pions = new TClonesArray("GenParticle");
    TClonesArray* reconstructed_neutral_pions = new TClonesArray("GenParticle");

    tree.Branch("BeamParticle", &beam_particles, 32000, 0);
    tree.Branch("GenTau", &generator_taus, 32000, 0);
    tree.Branch("GenChargedPion", &generator_charged_pions, 32000, 0);
    tree.Branch("GenNeutralPion", &generator_neutral_pions, 32000, 0);
    tree.Branch("ChargedPion", &reconstructed_charged_pions, 32000, 0);
    tree.Branch("NeutralPion", &reconstructed_neutral_pions, 32000, 0);

    constexpr double beam_energy = 2.13;
    constexpr double tau_momentum = 1.1745929254001146;
    constexpr double pion_boost_momentum = 0.5909200401438414;
    constexpr double pion_rest_momentum = 0.8829484500284907;
    constexpr double pion_energy = 1.0715709743251098;

    append_fixture_particle(*beam_particles, 0, -11, 1, 0.0, 0.0, beam_energy, beam_energy);
    append_fixture_particle(*beam_particles, 1, 11, -1, 0.0, 0.0, -beam_energy, beam_energy);
    append_fixture_particle(*generator_taus, 0, -15, 1, -tau_momentum, 0.0, 0.0, beam_energy);
    append_fixture_particle(*generator_taus, 1, 15, -1, tau_momentum, 0.0, 0.0, beam_energy);

    append_fixture_particle(*generator_charged_pions, 0, 211, 1, -pion_boost_momentum, 0.0,
                            pion_rest_momentum, pion_energy);
    append_fixture_particle(*generator_charged_pions, 1, -211, -1, pion_boost_momentum, pion_rest_momentum, 0.0,
                            pion_energy);
    append_fixture_particle(*reconstructed_charged_pions, 0, 211, 1, -pion_boost_momentum, 0.0,
                            pion_rest_momentum, pion_energy);
    append_fixture_particle(*reconstructed_charged_pions, 1, -211, -1, pion_boost_momentum, pion_rest_momentum, 0.0,
                            pion_energy);

    tree.Fill();
    tree.Write();
    file->Write();
    file->Close();

    delete beam_particles;
    delete generator_taus;
    delete generator_charged_pions;
    delete generator_neutral_pions;
    delete reconstructed_charged_pions;
    delete reconstructed_neutral_pions;
}

}  // namespace tauamp::delphes::test

#endif  // TAUAMP_TESTS_DELPHES_PION_PAIR_DELPHES_FIXTURE_H_
