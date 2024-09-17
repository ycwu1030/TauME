#include "LHEF.h"
#include "tauoo.h"

int main(int argc, char const *argv[]) {
    TChain *ch = new TChain("LHEF");
    for (size_t i = 2; i < argc; i++) {
        ch->Add(argv[i]);
    }

    LHEF *lf = new LHEF(ch);

    ME_EDM_Re_t<TauDecay_pi<false>, TauDecay_pi<true>> ME;

    dlv_t p_ep;
    dlv_t p_em;
    dlv_t p_taup;
    dlv_t p_taum;
    dlv_t p_nutau_p;
    dlv_t p_nutau_m;
    dlv_t p_pip;
    dlv_t p_pim;

    TFile *ff = new TFile(argv[1], "RECREATE");
    TTree *t1 = new TTree("TAUME", "TAUME");

    double ME0;
    double ME1;
    double ME2;
    double OO;
    t1->Branch("ME0", &ME0);
    t1->Branch("ME1", &ME1);
    t1->Branch("ME2", &ME2);
    t1->Branch("OO", &OO);
    int nevents = ch->GetEntries();
    for (size_t i = 0; i < nevents; i++) {
        ch->GetEntry(i);
        for (size_t ip = 0; ip < lf->Particle_size; ip++) {
            int pid = lf->Particle_PID[ip];
            double E = lf->Particle_E[ip];
            double px = lf->Particle_Px[ip];
            double py = lf->Particle_Py[ip];
            double pz = lf->Particle_Pz[ip];
            if (pid == 11) {
                p_em.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == -11) {
                p_ep.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == 15) {
                p_taum.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == -15) {
                p_taup.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == 16) {
                p_nutau_m.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == -16) {
                p_nutau_p.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == 211) {
                p_pip.SetPxPyPzE(px, py, pz, E);
            }
            if (pid == -211) {
                p_pim.SetPxPyPzE(px, py, pz, E);
            }
        }
        ME.set_taum_decay_momenta({p_nutau_m, p_pim});
        ME.set_taup_decay_momenta({p_nutau_p, p_pip});
        ME.set_production_momenta({p_ep, p_em, p_taup, p_taum});
        ME0 = ME.get_ME2_SM();
        ME1 = ME.get_ME2_Interference();
        ME2 = ME.get_ME2_BSM();
        OO = ME1 / ME0;
        t1->Fill();
    }
    t1->Write();
    ff->Close();
    return 0;
}
