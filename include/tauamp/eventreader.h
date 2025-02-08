#ifndef TAU_OO_IO_EVENT_READER_H_
#define TAU_OO_IO_EVENT_READER_H_

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include "classes/DelphesClasses.h"
#include "tauamp/constants.h"
#include "tauamp/utilities.h"

#ifndef DEBUG
#define DEBUG 0
#endif
namespace tauamp {

typedef enum { Pion = 0, Rho = 1, A1_3prong = 2, A1_1prong = 3 } TauDecayMode_t;

template <TauDecayMode_t Mode>
constexpr int get_n_charged_pion() {
    return 1;
}

template <>
constexpr int get_n_charged_pion<A1_3prong>() {
    return 3;
}

template <TauDecayMode_t Mode>
constexpr int get_n_neutral_pion() {
    return 0;
}
template <>
constexpr int get_n_neutral_pion<Rho>() {
    return 1;
}
template <>
constexpr int get_n_neutral_pion<A1_1prong>() {
    return 2;
}

typedef enum { PARTON = 0, DETECTOR = 1 } ReadingLevel_t;

template <TauDecayMode_t Mode_m = Rho, TauDecayMode_t Mode_p = Rho, ReadingLevel_t level = DETECTOR>
class EventReader_t {
public:
    static constexpr int k_Beam_PID = 11;
    static constexpr int k_N_ChargedPion_TauM = get_n_charged_pion<Mode_m>();
    static constexpr int k_N_ChargedPion_TauP = get_n_charged_pion<Mode_p>();
    static constexpr int k_N_NeutralPion_TauM = get_n_neutral_pion<Mode_m>();
    static constexpr int k_N_NeutralPion_TauP = get_n_neutral_pion<Mode_p>();
    static constexpr int k_N_Pion_TauM = k_N_ChargedPion_TauM + k_N_NeutralPion_TauM;
    static constexpr int k_N_Pion_TauP = k_N_ChargedPion_TauP + k_N_NeutralPion_TauP;
    static constexpr ReadingLevel_t k_Read_Level = level;
    typedef struct {
        std::vector<dlv_t> p_production;
        std::vector<dlv_t> p_tau_m_decay;
        std::vector<dlv_t> p_tau_p_decay;
    } MomentaList_t;
    // * production orders:
    // *     e+ e- tau+ tau-
    // * decay products orders:
    // *     nu pions
    // * pion orders:
    // * For Pion: pi0
    // * For Rho: pic pi0
    // * For A1_3prong:
    // *     for tau-: pi- pi- pi+
    // *     for tau+: pi+ pi+ pi-
    // * For A1_1prong: pic pi0 pi0

    EventReader_t(TChain *c) : m_reader(c), m_p_pion_tau_m(k_N_Pion_TauM), m_p_pion_tau_p(k_N_Pion_TauP) {
        m_branch_GenTau = m_reader.UseBranch("GenTau");
        m_branch_GenChargedPion = m_reader.UseBranch("GenChargedPion");
        m_branch_GenNeutralPion = m_reader.UseBranch("GenNeutralPion");
        m_branch_GenNeutrino = m_reader.UseBranch("GenNeutrino");
        m_branch_ChargedPion = m_reader.UseBranch("ChargedPion");
        m_branch_NeutralPion = m_reader.UseBranch("NeutralPion");
        m_branch_BeamParticle = m_reader.UseBranch("BeamParticle");
        m_rnd = new TRandom3();
    }
    ~EventReader_t() { delete m_rnd; }

    int get_entries() const { return m_reader.GetEntries(); };
    bool read_entry(int i) {
        bool good = m_reader.ReadEntry(i);
        if (!good) return good;
        return m_read_event();
        // return m_analysis_event();
    }

    MomentaList_t get_momenta_random() const {
        int tmp = m_rnd->Integer(2);
        return get_momenta(tmp);
    }
    MomentaList_t get_momenta(int i) const {
        std::vector<dlv_t> production(4);
        std::vector<dlv_t> decays_p(1 + k_N_Pion_TauP);
        std::vector<dlv_t> decays_m(1 + k_N_Pion_TauM);
        production[0] = m_p_beam_p;
        production[1] = m_p_beam_m;
        switch (i) {
            default:
            case 0:
                production[2] = m_p_tau_p_sol_1;
                production[3] = m_p_tau_m_sol_1;
                decays_p[0] = m_p_nu_tau_p_sol_1;
                decays_m[0] = m_p_nu_tau_m_sol_1;
                break;
            case 1:
                production[2] = m_p_tau_p_sol_2;
                production[3] = m_p_tau_m_sol_2;
                decays_p[0] = m_p_nu_tau_p_sol_2;
                decays_m[0] = m_p_nu_tau_m_sol_2;
                break;
        }
        std::copy(m_p_pion_tau_p.begin(), m_p_pion_tau_p.end(), decays_p.begin() + 1);
        std::copy(m_p_pion_tau_m.begin(), m_p_pion_tau_m.end(), decays_m.begin() + 1);
        return {production, decays_m, decays_p};
    }

private:
    ExRootTreeReader m_reader;
    TClonesArray *m_branch_GenTau;
    TClonesArray *m_branch_GenChargedPion;
    TClonesArray *m_branch_GenNeutralPion;
    TClonesArray *m_branch_GenNeutrino;
    TClonesArray *m_branch_ChargedPion;
    TClonesArray *m_branch_NeutralPion;
    TClonesArray *m_branch_BeamParticle;

    TRandom *m_rnd;

    // bool m_analysis_event() {
    //     if constexpr (k_Read_Level == DETECTOR) {  // C++17
    //         return m_read_detector_level_event();
    //     } else {
    //         return m_read_parton_level_event();
    //     }
    // }

    bool m_read_beam_particle() {
        int n_beam_particle = m_branch_BeamParticle->GetEntriesFast();
        if (n_beam_particle != 2) return false;
        // * Get Beam Particles
        GenParticle *beam1 = (GenParticle *)m_branch_BeamParticle->At(0);
        GenParticle *beam2 = (GenParticle *)m_branch_BeamParticle->At(1);
        int beamc1 = beam1->Charge;
        int beamc2 = beam2->Charge;
        if (beamc1 * beamc2 >= 0) return false;
        GenParticle *beam_m = beamc1 < 0 ? beam1 : beam2;
        GenParticle *beam_p = beamc1 > 0 ? beam1 : beam2;

        m_p_beam_m.SetPxPyPzE(beam_m->Px, beam_m->Py, beam_m->Pz, beam_m->E);
        m_p_beam_p.SetPxPyPzE(beam_p->Px, beam_p->Py, beam_p->Pz, beam_p->E);
        return true;
    }

    bool m_read_charged_pion(int &tau_m_ID, int &tau_p_ID) {
        // * Only for the case where there is only one charged pion in the decay products
        TClonesArray *branch;
        if constexpr (k_Read_Level == DETECTOR) {
            branch = m_branch_ChargedPion;
        } else {
            branch = m_branch_GenChargedPion;
        }
        int n_charged_pion = branch->GetEntriesFast();
        if (n_charged_pion != 2) return false;

        GenParticle *pic1 = (GenParticle *)branch->At(0);
        GenParticle *pic2 = (GenParticle *)branch->At(1);
        int charge1 = pic1->Charge;
        int charge2 = pic2->Charge;
        if (charge1 * charge2 >= 0) return false;
        GenParticle *pip = charge1 > 0 ? pic1 : pic2;
        GenParticle *pim = charge1 < 0 ? pic1 : pic2;
        tau_m_ID = pim->M1;
        tau_p_ID = pip->M1;

        m_p_pion_tau_m[0].SetPxPyPzE(pim->Px, pim->Py, pim->Pz, pim->E);
        m_p_pion_tau_p[0].SetPxPyPzE(pip->Px, pip->Py, pip->Pz, pip->E);
        return true;
    }

    bool m_read_tau() {}
    bool m_read_neutrino(int tau_m_id, int tau_p_id) {
        // * Reading truth neutrinos
        GenParticle *nu1 = (GenParticle *)m_branch_GenNeutrino->At(0);
        GenParticle *nu2 = (GenParticle *)m_branch_GenNeutrino->At(1);
        GenParticle *nup;
        GenParticle *num;
        int nu1_M1 = nu1->M1;
        int nu2_M1 = nu2->M1;
        if (nu1_M1 == tau_p_id && nu2_M1 == tau_m_id) {
            nup = nu1;
            num = nu2;
        } else if (nu1_M1 == tau_m_id && nu2_M1 == tau_p_id) {
            nup = nu2;
            num = nu1;
        } else {
            return false;
        }
        m_p_nu_tau_m_gen.SetPxPyPzE(num->Px, num->Py, num->Pz, num->E);
        m_p_nu_tau_p_gen.SetPxPyPzE(nup->Px, nup->Py, nup->Pz, nup->E);

        m_p_nu_tau_m_sol_1 = m_p_nu_tau_m_gen;
        m_p_nu_tau_p_sol_1 = m_p_nu_tau_p_gen;
        m_p_nu_tau_m_sol_2 = m_p_nu_tau_m_gen;
        m_p_nu_tau_p_sol_2 = m_p_nu_tau_p_gen;

        return true;
    }
    bool reconstruct_neutrinos(double Ecm, TLorentzVector &p_h_m, TLorentzVector &p_h_p) {
        double Etau = Ecm / 2.0;
        double ptau = Ecm / 2.0 * sqrt(1.0 - 4.0 * MTAU * MTAU / Ecm / Ecm);
        TLorentzVector pcm(0, 0, 0, Ecm);

        // TLorentzVector p_h_m = pi0m->P4() + pim->P4();
        // TLorentzVector p_h_p = pi0p->P4() + pip->P4();
        TLorentzVector p_miss = pcm - (p_h_m + p_h_p);
        double pm = p_miss.P();
        double thetam = p_miss.Theta();
        double phim = p_miss.Phi();

        double Enu = Ecm / 2.0 - p_h_m.E();
        double Enubar = Ecm / 2.0 - p_h_p.E();
        if (Enu < 0 || Enubar < 0) return false;
        if constexpr (DEBUG) {
            std::cout << "p_miss: (" << p_miss.Px() << "," << p_miss.Py() << "," << p_miss.Pz() << "," << p_miss.E()
                      << ")" << std::endl;
        }

        // * Go to p_miss frame: where p_miss in in the z direction
        // * By first rotation around z-axis by -phim angle
        // * then rotation around y-axis by -thetam angle
        p_h_m.RotateZ(-phim);
        p_h_m.RotateY(-thetam);
        p_h_p.RotateZ(-phim);
        p_h_p.RotateY(-thetam);
        p_miss.RotateZ(-phim);
        p_miss.RotateY(-thetam);

        if constexpr (DEBUG) {
            std::cout << "p_miss mag: " << pm << std::endl;
            std::cout << "p_miss: " << p_miss.E() << "," << p_miss.Px() << "," << p_miss.Py() << "," << p_miss.Pz()
                      << std::endl;
            std::cout << "p_h_m: " << p_h_m.E() << "," << p_h_m.Px() << "," << p_h_m.Py() << "," << p_h_m.Pz()
                      << std::endl;
            std::cout << "p_h_p: " << p_h_p.E() << "," << p_h_p.Px() << "," << p_h_p.Py() << "," << p_h_p.Pz()
                      << std::endl;
        }

        double cth = (pm * pm + Enu * Enu - Enubar * Enubar) / (2.0 * pm * Enu);
        if (abs(cth) > 1.0) return false;
        double theta = acos(cth);
        double sth = sin(theta);

        double z = Enu * cth;
        double r0 = Enu * sth;
        double x0 = p_h_m.Px();
        double y0 = p_h_m.Py();
        double z02 = ptau * ptau - pow(p_h_m.Pz() + z, 2);

        double ax = 4.0 * (x0 * x0 + y0 * y0);
        double bx = 4.0 * x0 * (r0 * r0 + x0 * x0 + y0 * y0 - z02);
        double cx = r0 * r0 * r0 * r0 + 2.0 * r0 * r0 * (x0 * x0 - y0 * y0 - z02) + pow(x0 * x0 + y0 * y0 - z02, 2);
        double Delta = bx * bx - 4.0 * ax * cx;
        if (Delta < 0) return false;
        double x_sol_1 = (-bx + sqrt(Delta)) / 2.0 / ax;
        double x_sol_2 = (-bx - sqrt(Delta)) / 2.0 / ax;

        double y_sol_1 = (z02 - r0 * r0 - y0 * y0 - x0 * x0 - 2.0 * x0 * x_sol_1) / 2.0 / y0;
        double y_sol_2 = (z02 - r0 * r0 - y0 * y0 - x0 * x0 - 2.0 * x0 * x_sol_2) / 2.0 / y0;

        // * Obtain the solution for neutrinos in the p_miss frame
        TLorentzVector p_nu_sol_1;
        p_nu_sol_1.SetXYZM(x_sol_1, y_sol_1, z, 0.0);
        TLorentzVector p_nubar_sol_1 = p_miss - p_nu_sol_1;

        TLorentzVector p_nu_sol_2;
        p_nu_sol_2.SetXYZM(x_sol_2, y_sol_2, z, 0.0);
        TLorentzVector p_nubar_sol_2 = p_miss - p_nu_sol_2;

        if constexpr (DEBUG) {
            std::cout << "p_nu_m_1: " << p_nu_sol_1.E() << "," << p_nu_sol_1.Px() << "," << p_nu_sol_1.Py() << ","
                      << p_nu_sol_1.Pz() << std::endl;
            std::cout << "p_nu_p_1: " << p_nubar_sol_1.E() << "," << p_nubar_sol_1.Px() << "," << p_nubar_sol_1.Py()
                      << "," << p_nubar_sol_1.Pz() << std::endl;

            std::cout << "p_nu_m_2: " << p_nu_sol_2.E() << "," << p_nu_sol_2.Px() << "," << p_nu_sol_2.Py() << ","
                      << p_nu_sol_2.Pz() << std::endl;
            std::cout << "p_nu_p_2: " << p_nubar_sol_2.E() << "," << p_nubar_sol_2.Px() << "," << p_nubar_sol_2.Py()
                      << "," << p_nubar_sol_2.Pz() << std::endl;
        }

        // * Rotate back to lab frame
        p_h_m.RotateY(thetam);
        p_h_m.RotateZ(phim);
        p_h_p.RotateY(thetam);
        p_h_p.RotateZ(phim);

        p_nu_sol_1.RotateY(thetam);
        p_nu_sol_1.RotateZ(phim);

        p_nubar_sol_1.RotateY(thetam);
        p_nubar_sol_1.RotateZ(phim);

        p_nu_sol_2.RotateY(thetam);
        p_nu_sol_2.RotateZ(phim);

        p_nubar_sol_2.RotateY(thetam);
        p_nubar_sol_2.RotateZ(phim);

        if constexpr (DEBUG) {
            p_miss.RotateY(thetam);
            p_miss.RotateZ(phim);
            std::cout << "p_miss: (" << p_miss.Px() << "," << p_miss.Py() << "," << p_miss.Pz() << "," << p_miss.E()
                      << ")" << std::endl;
        }

        // * Storing Neutrinos:
        m_p_nu_tau_m_sol_1.SetPxPyPzE(p_nu_sol_1.Px(), p_nu_sol_1.Py(), p_nu_sol_1.Pz(), p_nu_sol_1.E());
        m_p_nu_tau_p_sol_1.SetPxPyPzE(p_nubar_sol_1.Px(), p_nubar_sol_1.Py(), p_nubar_sol_1.Pz(), p_nubar_sol_1.E());

        m_p_nu_tau_m_sol_2.SetPxPyPzE(p_nu_sol_2.Px(), p_nu_sol_2.Py(), p_nu_sol_2.Pz(), p_nu_sol_2.E());
        m_p_nu_tau_p_sol_2.SetPxPyPzE(p_nubar_sol_2.Px(), p_nubar_sol_2.Py(), p_nubar_sol_2.Pz(), p_nubar_sol_2.E());

        return true;
    }

    bool m_read_event();

    // bool m_read_detector_level_event() { return m_read_detector_event<Mode_m, Mode_p>(); }

    dlv_t m_p_beam_m;
    dlv_t m_p_beam_p;
    dlv_t m_p_tau_m_sol_1;
    dlv_t m_p_tau_p_sol_1;
    dlv_t m_p_nu_tau_m_sol_1;
    dlv_t m_p_nu_tau_p_sol_1;

    dlv_t m_p_tau_m_sol_2;
    dlv_t m_p_tau_p_sol_2;
    dlv_t m_p_nu_tau_m_sol_2;
    dlv_t m_p_nu_tau_p_sol_2;

    dlv_t m_p_nu_tau_m_gen;
    dlv_t m_p_nu_tau_p_gen;

    std::vector<dlv_t> m_p_pion_tau_m;
    std::vector<dlv_t> m_p_pion_tau_p;
};

template <>
bool EventReader_t<Pion, Pion, PARTON>::m_read_event() {
    // * Get Beam Particles
    bool good = m_read_beam_particle();
    if (!good) return false;

    // * For pi+ pi-
    int tau_m_ID, tau_p_ID;
    good = m_read_charged_pion(tau_m_ID, tau_p_ID);
    if (!good) return false;

    // * For neutrino
    good = m_read_neutrino(tau_m_ID, tau_p_ID);  // * Read true neutrino
    if (!good) return false;

    if constexpr (DEBUG) {
        // * Checking whether the calculation of the neutrino can match the true neutrino
        double Ecm = m_p_beam_m.E() + m_p_beam_p.E();
        TLorentzVector p_h_m;
        TLorentzVector p_h_p;

        p_h_m.SetPxPyPzE(m_p_pion_tau_m[0].Px(), m_p_pion_tau_m[0].Py(), m_p_pion_tau_m[0].Pz(), m_p_pion_tau_m[0].E());
        p_h_p.SetPxPyPzE(m_p_pion_tau_p[0].Px(), m_p_pion_tau_p[0].Py(), m_p_pion_tau_p[0].Pz(), m_p_pion_tau_p[0].E());

        good = reconstruct_neutrinos(Ecm, p_h_m, p_h_p);
        if (!good) return false;
        std::cout << "True nu: " << m_p_nu_tau_m_gen << "," << "Reco nu: " << m_p_nu_tau_m_sol_1 << ","
                  << m_p_nu_tau_m_sol_2 << std::endl;
        std::cout << "True nubar: " << m_p_nu_tau_p_gen << "," << "Reco nubar: " << m_p_nu_tau_p_sol_1 << ","
                  << m_p_nu_tau_p_sol_2 << std::endl;
    }

    m_p_tau_m_sol_1 = m_p_nu_tau_m_sol_1 + m_p_pion_tau_m[0];
    m_p_tau_p_sol_1 = m_p_nu_tau_p_sol_1 + m_p_pion_tau_p[0];

    m_p_tau_m_sol_2 = m_p_nu_tau_m_sol_2 + m_p_pion_tau_m[0];
    m_p_tau_p_sol_2 = m_p_nu_tau_p_sol_2 + m_p_pion_tau_p[0];

    return true;
}

template <>
bool EventReader_t<Rho, Rho, PARTON>::m_read_event() {
    // * Get Beam Particles
    bool good = m_read_beam_particle();
    if (!good) return false;

    // * For Rho, Rho
    // ** Find charged pion for each tau, using charge to identify the origin
    int tau_p_id, tau_m_id;
    good = m_read_charged_pion(tau_m_id, tau_p_id);
    if (!good) {
        std::cout << "Charged Pion Error" << std::endl;
        return false;
    }

    // ** The neutrinos
    good = m_read_neutrino(tau_m_id, tau_p_id);
    if (!good) {
        std::cout << "Neutrino Error" << std::endl;
        return false;
    }

    // ** Neutral Pion, using the mother ID
    int n_neutral_pion = m_branch_GenNeutralPion->GetEntriesFast();
    if (n_neutral_pion != 2) {
        std::cout << "Neutral Pion Number Error" << std::endl;
        return false;
    }
    GenParticle *pi01 = (GenParticle *)m_branch_GenNeutralPion->At(0);
    GenParticle *pi02 = (GenParticle *)m_branch_GenNeutralPion->At(1);
    GenParticle *pi0p;
    GenParticle *pi0m;
    int pi01_M1 = pi01->M1;
    int pi02_M1 = pi02->M1;
    if (pi01_M1 == tau_p_id && pi02_M1 == tau_m_id) {
        pi0p = pi01;
        pi0m = pi02;
    } else if (pi01_M1 == tau_m_id && pi02_M1 == tau_p_id) {
        pi0p = pi02;
        pi0m = pi01;
    } else {
        std::cout << "Neutral Pion Matching Error" << std::endl;
        return false;
    }
    m_p_pion_tau_m[1].SetPxPyPzE(pi0m->Px, pi0m->Pz, pi0m->Pz, pi0m->E);
    m_p_pion_tau_p[1].SetPxPyPzE(pi0p->Px, pi0p->Py, pi0p->Pz, pi0p->E);

    // * Storing all the momenta

    m_p_tau_m_sol_1 = m_p_nu_tau_m_sol_1 + m_p_pion_tau_m[0] + m_p_pion_tau_m[1];
    m_p_tau_p_sol_1 = m_p_nu_tau_p_sol_1 + m_p_pion_tau_p[0] + m_p_pion_tau_p[1];

    m_p_tau_m_sol_2 = m_p_nu_tau_m_sol_2 + m_p_pion_tau_m[0] + m_p_pion_tau_m[1];
    m_p_tau_p_sol_2 = m_p_nu_tau_p_sol_2 + m_p_pion_tau_p[0] + m_p_pion_tau_p[1];

    return true;
}

template <>
bool EventReader_t<Pion, Pion, DETECTOR>::m_read_event() {
    bool good = m_read_beam_particle();
    if (!good) return false;
    double Ecm = m_p_beam_m.E() + m_p_beam_p.E();

    int tau_m_ID, tau_p_ID;
    good = m_read_charged_pion(tau_m_ID, tau_p_ID);
    if (!good) return false;

    // * Calculating the neutrino momentum
    TLorentzVector p_h_m;
    TLorentzVector p_h_p;

    p_h_m.SetPxPyPzE(m_p_pion_tau_m[0].Px(), m_p_pion_tau_m[0].Py(), m_p_pion_tau_m[0].Pz(), m_p_pion_tau_m[0].E());
    p_h_p.SetPxPyPzE(m_p_pion_tau_p[0].Px(), m_p_pion_tau_p[0].Py(), m_p_pion_tau_p[0].Pz(), m_p_pion_tau_p[0].E());
    good = reconstruct_neutrinos(Ecm, p_h_m, p_h_p);
    if (!good) return false;

    m_p_tau_m_sol_1 = m_p_nu_tau_m_sol_1 + m_p_pion_tau_m[0];
    m_p_tau_p_sol_1 = m_p_nu_tau_p_sol_1 + m_p_pion_tau_p[0];

    m_p_tau_m_sol_2 = m_p_nu_tau_m_sol_2 + m_p_pion_tau_m[0];
    m_p_tau_p_sol_2 = m_p_nu_tau_p_sol_2 + m_p_pion_tau_p[0];

    return true;
}

template <>
bool EventReader_t<Rho, Rho, DETECTOR>::m_read_event() {
    bool good = m_read_beam_particle();
    if (!good) return false;
    double Ecm = m_p_beam_m.E() + m_p_beam_p.E();

    int tau_m_ID, tau_p_ID;
    good = m_read_charged_pion(tau_m_ID, tau_p_ID);
    if (!good) return false;

    // * Pairing the pi0s
    // * two possible pairings:
    // * 1. pip + pi01 and pim + pi02
    // * 2. pip + pi02 and pim + pi01
    // * We use two information to determine the pairing
    // * invariant mass and dR
    // * Maybe more rely on invariant mass but not dR
    // * I will mainly determine the pairing by the invariant mass,
    // * only when the invariant mass gives the similar results I will use dR information
    TLorentzVector v_pip;
    v_pip.SetPxPyPzE(m_p_pion_tau_p[0].Px(), m_p_pion_tau_p[0].Py(), m_p_pion_tau_p[0].Pz(), m_p_pion_tau_p[0].E());
    TLorentzVector v_pim;
    v_pim.SetPxPyPzE(m_p_pion_tau_m[0].Px(), m_p_pion_tau_m[0].Py(), m_p_pion_tau_m[0].Pz(), m_p_pion_tau_m[0].E());

    int n_neutral_pion = m_branch_NeutralPion->GetEntriesFast();
    if (n_neutral_pion != k_N_NeutralPion_TauM + k_N_NeutralPion_TauP) return false;
    GenParticle *pi01 = (GenParticle *)m_branch_NeutralPion->At(0);
    GenParticle *pi02 = (GenParticle *)m_branch_NeutralPion->At(1);
    TLorentzVector v_pi01 = pi01->P4();
    TLorentzVector v_pi02 = pi02->P4();
    double minv11 = (v_pip + v_pi01).M();
    double minv12 = (v_pim + v_pi02).M();

    double minv21 = (v_pip + v_pi02).M();
    double minv22 = (v_pim + v_pi01).M();

    double chisq1 = pow((minv11 - MRHO) / GAMMARHO, 2) + pow((minv12 - MRHO) / GAMMARHO, 2);
    double chisq2 = pow((minv21 - MRHO) / GAMMARHO, 2) + pow((minv22 - MRHO) / GAMMARHO, 2);
    bool allclose = (chisq1 <= 2.0 && chisq2 <= 2.0);

    if constexpr (DEBUG) {
        std::cout << "Minv 1: " << minv11 << " " << minv12 << std::endl;
        std::cout << "Minv 2: " << minv21 << " " << minv22 << std::endl;
        std::cout << "chisq1: " << chisq1 << " chisq2: " << chisq2 << std::endl;
    }

    int good_choice = 1;  // *  1 or 2
    if (allclose) {
        // * In either case, we have good invariant mass, then we resort to dR
        double dR11 = v_pip.DeltaR(v_pi01);
        double dR12 = v_pim.DeltaR(v_pi02);
        double dR1 = dR11 + dR12;

        double dR21 = v_pip.DeltaR(v_pi02);
        double dR22 = v_pim.DeltaR(v_pi01);
        double dR2 = dR21 + dR22;

        good_choice = dR1 < dR2 ? 1 : 2;
    } else {
        good_choice = chisq1 < chisq2 ? 1 : 2;
    }
    GenParticle *pi0p = good_choice == 1 ? pi01 : pi02;
    GenParticle *pi0m = good_choice == 1 ? pi02 : pi01;
    m_p_pion_tau_m[1].SetPxPyPzE(pi0m->Px, pi0m->Py, pi0m->Pz, pi0m->E);
    m_p_pion_tau_p[1].SetPxPyPzE(pi0p->Px, pi0p->Py, pi0p->Pz, pi0p->E);

    // * Calculating the neutrino momentum
    TLorentzVector p_h_m = pi0m->P4() + v_pim;
    TLorentzVector p_h_p = pi0p->P4() + v_pip;
    good = reconstruct_neutrinos(Ecm, p_h_m, p_h_p);
    if (!good) return false;

    if constexpr (DEBUG) {
        double p_nu_mag = m_p_nu_tau_m_gen.P();
        double p_nu_sol_1_mag = m_p_nu_tau_m_sol_1.P();
        double p_nu_sol_2_mag = m_p_nu_tau_m_sol_2.P();
        double dp1 = abs(p_nu_sol_1_mag - p_nu_mag);
        double dp2 = abs(p_nu_sol_2_mag - p_nu_mag);
        int good_p_mag = dp1 < dp2 ? 1 : 2;
        double dp_min = dp1 < dp2 ? dp1 : dp2;

        double dR1 = ROOT::Math::VectorUtil::DeltaR(m_p_nu_tau_m_gen, m_p_nu_tau_m_sol_1);
        double dR2 = ROOT::Math::VectorUtil::DeltaR(m_p_nu_tau_m_gen, m_p_nu_tau_m_sol_2);
        int good_dR = dR1 < dR2 ? 1 : 2;
        double dR_min = dR1 < dR2 ? dR1 : dR2;

        std::cout << "mag: " << good_p_mag << " dp_min = " << dp_min << std::endl;
        std::cout << "dR: " << good_dR << " dR_min = " << dR_min << std::endl;
    }

    // * Storing Tau:
    m_p_tau_m_sol_1 = m_p_nu_tau_m_sol_1 + m_p_pion_tau_m[0] + m_p_pion_tau_m[1];
    m_p_tau_p_sol_1 = m_p_nu_tau_p_sol_1 + m_p_pion_tau_p[0] + m_p_pion_tau_p[1];

    m_p_tau_m_sol_2 = m_p_nu_tau_m_sol_2 + m_p_pion_tau_m[0] + m_p_pion_tau_m[1];
    m_p_tau_p_sol_2 = m_p_nu_tau_p_sol_2 + m_p_pion_tau_p[0] + m_p_pion_tau_p[1];

    return true;
}

// class RhoReader_t {
// public:
//     RhoReader_t(TChain *c) : m_reader(c) {
//         m_branch_GenChargedPion = m_reader.UseBranch("GenChargedPion");
//         m_branch_GenNeutralPion = m_reader.UseBranch("GenNeutralPion");
//         m_branch_GenNeutrino = m_reader.UseBranch("GenNeutrino");
//         m_branch_ChargedPion = m_reader.UseBranch("ChargedPion");
//         m_branch_NeutralPion = m_reader.UseBranch("NeutralPion");
//         m_branch_BeamParticle = m_reader.UseBranch("BeamParticle");
//     }
//     ~RhoReader_t() {}

// private:
//     ExRootTreeReader m_reader;
//     TClonesArray *m_branch_GenChargedPion;
//     TClonesArray *m_branch_GenNeutralPion;
//     TClonesArray *m_branch_GenNeutrino;
//     TClonesArray *m_branch_ChargedPion;
//     TClonesArray *m_branch_NeutralPion;
//     TClonesArray *m_branch_BeamParticle;
// };

};  // namespace tauamp

#endif  // TAU_OO_IO_EVENT_READER_H_
