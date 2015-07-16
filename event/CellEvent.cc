#include "CellEvent.h"

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>

#include "TNamed.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

using namespace std;

CellEvent::CellEvent(){}

CellEvent::CellEvent(const char* filename)
{
    dbPDG = new TDatabasePDG();
    thresh_KE = 1; // MeV

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

    trackIndex.clear();
    trackParents.clear();
    trackChildren.clear();
    trackSiblings.clear();
}

void CellEvent::GetEntry(int entry)
{
    Reset();
    simTree->GetEntry(entry);
    ProcessTracks();
    // reco_trackPosition->Print();
    // TClonesArray *pos = (TClonesArray*)(*reco_trackPosition)[reco_nTrack-1];
    // pos->Print();

}


//----------------------------------------------------------------
void CellEvent::ProcessTracks()
{
    // map track id to track index in the array
    for (int i=0; i<mc_Ntrack; i++) {
        trackIndex[mc_id[i]] = i;
    }

    // in trackParents, trackChildren, trackSiblings vectors, store track index (not track id)
    for (int i=0; i<mc_Ntrack; i++) {
        // currently, parent size == 1;
        // for primary particle, parent id = 0;
        vector<int> parents;
        if ( !IsPrimary(i) ) {
            parents.push_back(trackIndex[mc_mother[i]]);
        }
        trackParents.push_back(parents); // primary track will have 0 parents

        vector<int> children;
        int nChildren = (*mc_daughters).at(i).size();
        for (int j=0; j<nChildren; j++) {
            children.push_back(trackIndex[(*mc_daughters).at(i).at(j)]);
        }
        trackChildren.push_back(children);

    }

    // siblings
    for (int i=0; i<mc_Ntrack; i++) {
        vector<int> siblings;
        if ( IsPrimary(i) ) {
            for (int j=0; j<mc_Ntrack; j++) {
                if( IsPrimary(j) ) {
                    siblings.push_back(j);
                }
            }
        }
        else {
            // siblings are simply children of the mother
            int mother = trackIndex[mc_mother[i]];
            int nSiblings = trackChildren.at(mother).size();
            for (int j=0; j<nSiblings; j++) {
                siblings.push_back(trackChildren.at(mother).at(j));
            }
        }
        trackSiblings.push_back(siblings);
    }

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

void CellEvent::PrintMCInfo()
{
    cout << "number of tracks:" << mc_Ntrack << endl;
    for (int i=0; i<mc_Ntrack; i++) {
        if (mc_pdg[i] == 11 || mc_pdg[i] == -11 || mc_pdg[i] > 1e9 || mc_pdg[i] == 22) continue;
        cout << i <<  ": track id = " << mc_id[i] << ", pdg = " << mc_pdg[i] << endl;
        cout << "parents: "; PrintVector(trackParents[i]); cout << endl;
        cout << "--------" << endl;
    }
}

void CellEvent::PrintVector(vector<int>& v)
{
    int size = v.size();
    if (size==0) return;

    for (int i=0; i<size-1; i++) {
        cout << v[i] << ", ";
    }
    cout << v[size-1];
}

double CellEvent::KE(float* momentum)
{
    TLorentzVector particle(momentum);
    return particle.E()-particle.M();
}

bool CellEvent::KeepMC(int i)
{
    double e = KE(mc_startMomentum[i])*1000;
    double thresh_KE_em = 5.; // MeV
    double thresh_KE_np = 10; // MeV
    if (mc_pdg[i]==22 || mc_pdg[i]==11 || mc_pdg[i]==-11) {
        if (e>=thresh_KE_em) return true;
        else return false;
    }
    else if (mc_pdg[i]==2112 || mc_pdg[i]==2212 || mc_pdg[i]>1e9) {
        if (e>=thresh_KE_np) return true;
        else return false;
    }
    return true;
}

bool CellEvent::DumpJSON(int id, ostream& out)
{
    int i = trackIndex[id];
    if (!KeepMC(i)) return false;

    int e = KE(mc_startMomentum[i])*1000;

    int nDaughter = (*mc_daughters).at(i).size();
    vector<int> saved_daughters;
    for (int j=0; j<nDaughter; j++) {
        int daughter_id = (*mc_daughters).at(i).at(j);
        // int e_daughter = KE(mc_startMomentum[ trackIndex[daughter_id] ])*1000;
        // if (e_daughter >= thresh_KE) {
        if ( KeepMC(trackIndex[daughter_id]) ) {
            saved_daughters.push_back(daughter_id);
        }
    }

    out << fixed << setprecision(1);
    out << "{";

    out << "\"id\":" << id << ",";
    out << "\"text\":" << "\"" << PDGName(mc_pdg[i]) << "  " << e << " MeV\",";
    out << "\"data\":{";
    out << "\"start\":[" << mc_startXYZT[i][0] << ", " <<  mc_startXYZT[i][1] << ", " << mc_startXYZT[i][2] << "],";
    out << "\"end\":[" << mc_endXYZT[i][0] << ", " <<  mc_endXYZT[i][1] << ", " << mc_endXYZT[i][2] << "]";
    out << "},";
    out << "\"children\":[";
    int nSavedDaughter = saved_daughters.size();
    if (nSavedDaughter == 0) {
        out << "],";
        out << "\"icon\":" << "\"jstree-file\"";
        out << "}";
        return true;
    }
    else {
        for (int j=0; j<nSavedDaughter; j++) {
            DumpJSON(saved_daughters.at(j), out);
            if (j!=nSavedDaughter-1) {
                out << ",";
            }
        }
        out << "]";
        out << "}";
        return true;
    }
}

void CellEvent::DumpJSON(ostream& out)
{
    out << "[";
    vector<int> primaries;
    for (int i=0; i<mc_Ntrack; i++) {
        if (IsPrimary(i)) {
            // int e = KE(mc_startMomentum[i])*1000;
            // if (e<thresh_KE) continue;
            if (KeepMC(i)) {
                primaries.push_back(i);
            }
        }
    }
    int size = primaries.size();
    // cout << size << endl;
    for (int i=0; i<size; i++) {
        if (DumpJSON(mc_id[primaries[i]], out) && i!=size-1) {
            out << ", ";
        }
    }

    out << "]";
}

TString CellEvent::PDGName(int pdg)
{
    TParticlePDG *p = dbPDG->GetParticle(pdg);
    if (p == 0) {
        if (pdg>1e9) {
            int z = (pdg - 1e9) / 10000;
            int a = (pdg - 1e9 - z*1e4) / 10;
            TString name;
            if (z == 18) name = "Ar";

            else if (z == 17) name = "Cl";
            else if (z == 19) name = "Ca";
            else if (z == 16) name = "S";
            else if (z == 15) name = "P";
            else if (z == 14) name = "Si";
            else if (z == 1) name = "H";
            else if (z == 2) name = "He";

            else return pdg;
            return Form("%s-%i", name.Data(), a);
        }
        return pdg;
    }
    else {
        return p->GetName();
    }
}
