#ifndef CELLEVENT_H
#define CELLEVENT_H

#include <vector>
#include <map>
#include <iostream>
#include "TString.h"

class TClonesArray;
class TTree;
class TFile;
class TDatabasePDG;

class CellEvent {
public:
    TFile *rootFile;
    TTree *simTree;

    enum LIMITS {
        MAX_TRACKS = 30000,
    };

    int nEvent;

    // simTree Leafs
    int eventNo;
    int runNo;
    int subRunNo;

    int raw_nChannel;  // number of raw hit channels (including noise)
    std::vector<int> *raw_channelId;
    TClonesArray *raw_wf;

    int calib_nChannel;  // number of hit channels
    std::vector<int> *calib_channelId;
    TClonesArray *calib_wf;

    int simide_size;
    std::vector<int>   *simide_channelIdY;
    std::vector<int>   *simide_trackId;
    std::vector<int>   *simide_tdc;
    std::vector<float> *simide_x;
    std::vector<float> *simide_y;
    std::vector<float> *simide_z;
    std::vector<float> *simide_numElectrons;

    int mc_Ntrack;  // number of tracks in MC
    int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
    int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
    int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
    int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
    float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
    float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
    float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
    float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > *mc_daughters;  // daughters id of this track; vector

    int mc_isnu; // is neutrino interaction
    int mc_nGeniePrimaries; // number of Genie primaries
    int mc_nu_pdg; // pdg code of neutrino
    int mc_nu_ccnc; // cc or nc
    int mc_nu_mode; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
    int mc_nu_intType; // interaction type
    int mc_nu_target; // target interaction
    int mc_hitnuc; // hit nucleon
    int mc_hitquark; // hit quark
    double mc_nu_Q2; // Q^2
    double mc_nu_W; // W
    double mc_nu_X; // X
    double mc_nu_Y; // Y
    double mc_nu_Pt; // Pt
    double mc_nu_Theta; // angle relative to lepton
    float mc_nu_pos[4];  // interaction position of nu
    float mc_nu_mom[4];  // interaction momentum of nu

    // ----- derived ---
    std::map<int, int> trackIndex;
    std::vector<std::vector<int> > trackParents;
    std::vector<std::vector<int> > trackChildren;
    std::vector<std::vector<int> > trackSiblings;
    TDatabasePDG *dbPDG;

    //-------------------------------------
    CellEvent();
    CellEvent(const char* filename);
    virtual ~CellEvent();

    TTree* Tree() { return simTree; }
    void GetEntry(int i);
    void Reset();
    void PrintInfo(int level=0);  // print the current event(entry) info
    void PrintMCInfo();  // print mc tracks info
    void InitBranchAddress();

    bool IsPrimary(int i) { return mc_mother[i] == 0 ; }
    void ProcessTracks();
    void PrintVector(std::vector<int>& v);
    double KE(float* momentum);  // KE
    TString PDGName(int pdg);
    double thresh_KE;
    bool KeepMC(int i);

    bool DumpJSON(int id, ostream& out);
    void DumpJSON(ostream& out = std::cout);


};

#endif
