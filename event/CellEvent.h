#ifndef CELLEVENT_H
#define CELLEVENT_H

#include <vector>
class TClonesArray;
class TTree;
class TFile;

class CellEvent {
public:
    TFile *rootFile;
    TTree *simTree;

    int nEvent;

    // simTree Leafs
    int eventNo;
    int runNo;
    int subRunNo;

    int calib_nChannel;  // number of hit channels
    std::vector<int> *calib_channelId;
    TClonesArray *calib_wf;

    //-------------------------------------
    CellEvent();
    CellEvent(const char* filename);
    virtual ~CellEvent();

    TTree* Tree() { return simTree; }
    void GetEntry(int i);
    void Reset();
    void PrintInfo(int level=0);  // print the current event(entry) info
    void InitBranchAddress();

};

#endif