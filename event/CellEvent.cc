#include "CellEvent.h"

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TClonesArray.h"

using namespace std;

CellEvent::CellEvent(){}

CellEvent::CellEvent(const char* filename)
{
    calib_channelId = new vector<int>;
    calib_wf = new TClonesArray();

    rootFile = new TFile(filename);
    simTree = (TTree*)rootFile->Get("/Event/Sim");
    nEvent = simTree->GetEntries();

    InitBranchAddress();
}


CellEvent::~CellEvent()
{
    delete calib_channelId;
    delete calib_wf;

    rootFile->Close();
    delete rootFile;    
}


void CellEvent::InitBranchAddress()
{
    simTree->SetBranchAddress("eventNo" , &eventNo);
    simTree->SetBranchAddress("runNo"   , &runNo);
    simTree->SetBranchAddress("subRunNo", &subRunNo);

    simTree->SetBranchAddress("calib_nChannel", &calib_nChannel);
    simTree->SetBranchAddress("calib_channelId", &calib_channelId);
    simTree->SetBranchAddress("calib_wf", &calib_wf);
}

void CellEvent::Reset()
{
    (*calib_channelId).clear();
    calib_wf->Clear();
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
    if (level>1) {
        for (int i=0; i<calib_nChannel; i++) {
            cout << calib_channelId->at(i) << " ";
        }
        cout << endl;
    }
    cout << endl;
}
