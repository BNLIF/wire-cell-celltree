#include "CellEvent.h"

#include <iostream>
#include <vector>

#include "TNamed.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TClonesArray.h"

using namespace std;

CellEvent::CellEvent(){}

CellEvent::CellEvent(const char* filename)
{
    raw_channelId = new vector<int>;
    raw_wf = new TClonesArray();

    calib_channelId = new vector<int>;
    calib_wf = new TClonesArray();

    simide_channelIdY = new vector<int>;
    simide_trackId = new vector<int>;
    simide_tdc = new vector<int>;
    simide_x = new vector<float>;
    simide_y = new vector<float>;
    simide_z = new vector<float>;
    simide_numElectrons = new vector<float>;

    mc_daughters = new std::vector<std::vector<int> >;  // daughters id of this track; vector

    rootFile = new TFile(filename);
    simTree = (TTree*)rootFile->Get("/Event/Sim");
    nEvent = simTree->GetEntries();

    InitBranchAddress();
}


CellEvent::~CellEvent()
{
    delete raw_channelId;
    delete raw_wf;

    delete calib_channelId;
    delete calib_wf;
    delete mc_daughters;

    rootFile->Close();
    delete rootFile;
}


void CellEvent::InitBranchAddress()
{
    simTree->SetBranchAddress("eventNo" , &eventNo);
    simTree->SetBranchAddress("runNo"   , &runNo);
    simTree->SetBranchAddress("subRunNo", &subRunNo);

    simTree->SetBranchAddress("raw_nChannel", &raw_nChannel);
    simTree->SetBranchAddress("raw_channelId", &raw_channelId);
    simTree->SetBranchAddress("raw_wf", &raw_wf);

    simTree->SetBranchAddress("calib_nChannel", &calib_nChannel);
    simTree->SetBranchAddress("calib_channelId", &calib_channelId);
    simTree->SetBranchAddress("calib_wf", &calib_wf);

    simTree->SetBranchAddress("simide_size", &simide_size);
    simTree->SetBranchAddress("simide_channelIdY", &simide_channelIdY);
    simTree->SetBranchAddress("simide_trackId", &simide_trackId);
    simTree->SetBranchAddress("simide_tdc", &simide_tdc);
    simTree->SetBranchAddress("simide_x", &simide_x);
    simTree->SetBranchAddress("simide_y", &simide_y);
    simTree->SetBranchAddress("simide_z", &simide_z);
    simTree->SetBranchAddress("simide_numElectrons", &simide_numElectrons);

    simTree->SetBranchAddress("mc_Ntrack"       , &mc_Ntrack);
    simTree->SetBranchAddress("mc_id"           , &mc_id);
    simTree->SetBranchAddress("mc_pdg"          , &mc_pdg);
    simTree->SetBranchAddress("mc_process"      , &mc_process);
    simTree->SetBranchAddress("mc_mother"       , &mc_mother);
    simTree->SetBranchAddress("mc_daughters"    , &mc_daughters);
    simTree->SetBranchAddress("mc_startXYZT"    , &mc_startXYZT);
    simTree->SetBranchAddress("mc_endXYZT"      , &mc_endXYZT);
    simTree->SetBranchAddress("mc_startMomentum", &mc_startMomentum);
    simTree->SetBranchAddress("mc_endMomentum"  , &mc_endMomentum);

    simTree->SetBranchAddress("mc_isnu"  , &mc_isnu);
    simTree->SetBranchAddress("mc_nGeniePrimaries"  , &mc_nGeniePrimaries);
    simTree->SetBranchAddress("mc_nu_pdg"  , &mc_nu_pdg);
    simTree->SetBranchAddress("mc_nu_ccnc"  , &mc_nu_ccnc);
    simTree->SetBranchAddress("mc_nu_mode"  , &mc_nu_mode);
    simTree->SetBranchAddress("mc_nu_intType"  , &mc_nu_intType);
    simTree->SetBranchAddress("mc_nu_target"  , &mc_nu_target);
    simTree->SetBranchAddress("mc_hitnuc"  , &mc_hitnuc);
    simTree->SetBranchAddress("mc_hitquark"  , &mc_hitquark);
    simTree->SetBranchAddress("mc_nu_Q2"  , &mc_nu_Q2);
    simTree->SetBranchAddress("mc_nu_W"  , &mc_nu_W);
    simTree->SetBranchAddress("mc_nu_X"  , &mc_nu_X);
    simTree->SetBranchAddress("mc_nu_Y"  , &mc_nu_Y);
    simTree->SetBranchAddress("mc_nu_Pt"  , &mc_nu_Pt);
    simTree->SetBranchAddress("mc_nu_X"  , &mc_nu_X);
    simTree->SetBranchAddress("mc_nu_pos"  , &mc_nu_pos);
    simTree->SetBranchAddress("mc_nu_mom"  , &mc_nu_mom);

}

void CellEvent::Reset()
{
    (*raw_channelId).clear();
    raw_wf->Clear();

    (*calib_channelId).clear();
    calib_wf->Clear();

    (*mc_daughters).clear();

    simide_channelIdY->clear();
    simide_trackId->clear();
    simide_tdc->clear();
    simide_x->clear();
    simide_y->clear();
    simide_z->clear();
    simide_numElectrons->clear();
}

void CellEvent::GetEntry(int entry)
{
    Reset();
    simTree->GetEntry(entry);

    // reco_trackPosition->Print();
    // TClonesArray *pos = (TClonesArray*)(*reco_trackPosition)[reco_nTrack-1];
    // pos->Print();

}

void CellEvent::PrintInfo(int level)
{
    TNamed* version = (TNamed*)rootFile->Get("version");
    cout << "data version: " << version->GetTitle() << endl;

    cout << "run/subRun/event (total) : "
        << runNo << "/"
        << subRunNo << "/"
        << eventNo-1 << " ("
        << nEvent << ")"
        << endl;
    cout << "nChannel: " << calib_nChannel << endl;
    if (level>0) {
        cout << "first channel: " << endl;
        TH1F *h = (TH1F*)calib_wf->At(0);
        for (int i=0; i<9600; i++) {
            if (fabs(h->GetBinContent(i))>0.1) {
                cout << h->GetBinContent(i) << " ";
            }
        }
        cout << endl;
    }

    if (level > 0) {
        cout << "MC tracks:" << mc_Ntrack;
        for (int i=0; i<mc_Ntrack; i++) {
            cout << "\n              id: " << mc_id[i];
            cout << "\n             pdg: " << mc_pdg[i];
            cout << "\n          mother: " << mc_mother[i];
            cout << "\n      Ndaughters: " << (*mc_daughters).at(i).size();
            cout << "\n      start XYZT: (" << mc_startXYZT[i][0] << ", " << mc_startXYZT[i][1] << ", " << mc_startXYZT[i][2] << ", " << mc_startXYZT[i][3] << ")";
            cout << "\n        end XYZT: (" << mc_endXYZT[i][0] << ", " << mc_endXYZT[i][1] << ", " << mc_endXYZT[i][2] << ", " << mc_endXYZT[i][3] << ")";
            cout << "\n  start momentum: (" << mc_startMomentum[i][0] << ", " << mc_startMomentum[i][1] << ", " << mc_startMomentum[i][2] << ", " << mc_startMomentum[i][3] << ")";
            cout << "\n    end momentum: (" << mc_endMomentum[i][0] << ", " << mc_endMomentum[i][1] << ", " << mc_endMomentum[i][2] << ", " << mc_endMomentum[i][3] << ")";
            cout << endl;
        }
    }
    
    if (level>1) {
        for (int i=0; i<calib_nChannel; i++) {
            cout << calib_channelId->at(i) << " ";
        }
        cout << endl;
    }
    cout << endl;
}
