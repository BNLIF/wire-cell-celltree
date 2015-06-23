// Read data from MC raw files and convert it into ROOT tree
// Chao Zhang (chao@bnl.gov) 5/13/2014

#ifndef CTree_Module
#define CTree_Module

// LArSoft includes
#include "Utilities/DetectorProperties.h"
#include "Utilities/GeometryUtilities.h"
// #include "Utilities/LArProperties.h"

#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "AnalysisBase/ParticleID.h"
#include "AnalysisBase/Calorimetry.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindMany.h"
#include "art/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"


// ROOT includes.
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TString.h"

// C++ Includes
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

#define MAX_PLANE 3
#define MAX_CHANNEL 8254
#define MAX_TRACKS 3000
#define MAX_HITS 30000
#define MAX_TRACKER 10
#define MAX_TRACKER_TRACKS 300
#define MAX_TRACKER_HITS 30000

using namespace std;

namespace microboone {

class CTree : public art::EDAnalyzer {
public:

    explicit CTree(fhicl::ParameterSet const& pset);
    virtual ~CTree();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void initOutput();
    void saveChannelWireMap();
    void saveWireGeometry();
    void saveChannelWireGeometry();
    void printGeometry();
    void printEvent();

    void processMC(const art::Event& evt);
    void processRaw(const art::Event& evt);
    void processCalib(const art::Event& evt);
    void processHits(const art::Event& evt);
    void processTrack(const string& trackLabel, const art::Event& evt, int trackerIndex);
    void reset();

    void InitProcessMap();
    double TrackLength(const recob::Track& track);

private:
    std::map<std::string, int> processMap;

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    std::string fCalibLabel;
    std::string fHitsModuleLabel;
    std::vector<std::string> fTrackModuleLabel;

    std::string fOutFileName;
    bool fSaveChannelWireMap;

    art::ServiceHandle<geo::Geometry> fGeom;
    // // art::ServiceHandle<util::LArProperties> larp;

    TFile *fOutFile;
    TTree *fGeoTree;
    TTree *fEventTree;

    // Geometry Tree Leafs 
    float  fTPC_x;  // TPC length in x
    float  fTPC_y;  // TPC length in y
    float  fTPC_z;  // TPC length in z
    int    fNplanes;
    int    fPlane_type[MAX_PLANE];  // plane type: 0 == induction, 1 == collection
    int    fPlane_view[MAX_PLANE];  // wire orientation: 0 == U, 1 == V, 2 == X
    double fPlane_wirepitch[MAX_PLANE];  // wire pitch of each plane
    double fPlane_wireangle[MAX_PLANE];  // wire angle of each plane in radian
    int    fPlane_wires[MAX_PLANE];  // number of wires in each plane
    int    fNchannels;
    int    fNOpDet;


    // Event Tree Leafs
    int fEvent;
    int fRun;
    int fSubRun;

    int fRaw_Nhit;
    int fRaw_channelId[MAX_CHANNEL];  // hit channel id; size == raw_Nhit
    int fRaw_baseline[MAX_CHANNEL]; // hit channel baseline (noise level)
    int fRaw_charge[MAX_CHANNEL];  // hit channel charge (simple alg); size == raw_Nhit
    int fRaw_time[MAX_CHANNEL];  // hit channel time (simple alg); size == raw_Nhit
    std::vector<std::vector<int> > fRaw_wfADC;
    std::vector<std::vector<int> > fRaw_wfTDC;

    int fCalib_Nhit;
    int fCalib_channelId[MAX_CHANNEL];  // hit channel id; size == fCalib_Nhit
    // FIXEME:: cannot save e.g std::vector<std::vector<float> > in ttree
    std::vector<std::vector<int> > fCalib_wfADC;  
    std::vector<std::vector<int> > fCalib_wfTDC;


    int fMC_Ntrack;  // number of tracks in MC
    int fMC_id[MAX_TRACKS];  // track id; size == mc_Ntrack
    int fMC_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
    int fMC_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
    int fMC_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
    float fMC_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
    float fMC_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
    float fMC_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
    float fMC_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > fMC_daughters;  // daughters id of this track; vector

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

    int    no_hits;                  //number of hits
    int    hit_channel[MAX_HITS];    //channel ID
    int    hit_plane[MAX_HITS];      // plane #
    float  hit_peakT[MAX_HITS];      //peak time
    float  hit_charge[MAX_HITS];     //charge (area)

    int trk_nTrack[MAX_TRACKER]; // no. of tracks of each tracker.
    std::vector<int> trk_nHit[MAX_TRACKER];  // no. of hits of each track
    std::vector<double> trk_length[MAX_TRACKER];  // length of each track
    std::vector<std::vector<double> > trk_start_xyz[MAX_TRACKER];  // position of start vertex
    std::vector<std::vector<double> > trk_end_xyz[MAX_TRACKER];  // position of end vertex
    std::vector<std::vector<double> > trk_start_dxyz[MAX_TRACKER];  // direction of start vertex
    std::vector<std::vector<double> > trk_end_dxyz[MAX_TRACKER];  // direction of end vertex

    std::vector<std::vector<double> > trk_points_x[MAX_TRACKER];  // position of all points on trajectory
    std::vector<std::vector<double> > trk_points_y[MAX_TRACKER];  // position of all points on trajectory
    std::vector<std::vector<double> > trk_points_z[MAX_TRACKER];  // position of all points on trajectory

    std::vector<double> trk_calo_KE[MAX_TRACKER][3];  // KE of calorimetry from each plan of each trk
    std::vector<double> trk_calo_range[MAX_TRACKER][3];  // Range of calorimetry from each plan of each trk
    std::vector<int> trk_calo_nHit[MAX_TRACKER][3];  // hits of calorimetry from each plan of each trk
    std::vector<std::vector<double> > trk_calo_dedx[MAX_TRACKER][3];  // dedx of calorimetry from each plan of each trk
    std::vector<std::vector<double> > trk_calo_dqdx[MAX_TRACKER][3];  // dqdx of calorimetry from each plan of each trk
    std::vector<std::vector<double> > trk_calo_resRange[MAX_TRACKER][3];  // residual range of calorimetry from each plan of each trk

    // int no_trackers; // number of trackers
    // int trk_id[MAX_TRACKER][MAX_TRACKER_TRACKS];  // id of identified tracks; 
    // #define MAX_TRACKER 10
    // #define MAX_TRACKER_TRACKS 300
    // #define MAX_TRACKER_HITS 30000

}; // class CTree


//-----------------------------------------------------------------------
CTree::CTree(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
    InitProcessMap();
    initOutput();
}

//-----------------------------------------------------------------------
CTree::~CTree()
{
}

//-----------------------------------------------------------------------
void CTree::reconfigure(fhicl::ParameterSet const& p){
    fRawDigitLabel = p.get<std::string>("RawDigitLabel");
    fCalibLabel = p.get<std::string>("CalibLabel");
    fHitsModuleLabel = p.get<std::string>("HitsModuleLabel");
    fTrackModuleLabel = p.get< std::vector<std::string> >("TrackModuleLabel");
    fOutFileName = p.get<std::string>("outFile");
    fSaveChannelWireMap = p.get<bool>("saveChannelWireMap");
}

//-----------------------------------------------------------------------
void CTree::initOutput()
{
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");

    // init Detector Geometry TTree
    TDirectory* subDir = fOutFile->mkdir("Detector");
    subDir->cd();
    fGeoTree = new TTree("Geometry", "Detector Geometry");
    
    fGeoTree->Branch("TPC_x", &fTPC_x); // TPC length in x
    fGeoTree->Branch("TPC_y", &fTPC_y); // TPC length in y
    fGeoTree->Branch("TPC_z", &fTPC_z); // TPC length in z
    
    fGeoTree->Branch("Nplanes"     , &fNplanes);  // number of wire planes in each TPC, == 3
    fGeoTree->Branch("plane_type"  , &fPlane_type, "plane_type[Nplanes]/I"); // plane type: 0 == induction, 1 == collection
    fGeoTree->Branch("plane_view"  , &fPlane_view, "plane_view[Nplanes]/I"); // wire orientation: 0 == U, 1 == V, 2 == Z
    fGeoTree->Branch("plane_wirepitch"  , &fPlane_wirepitch, "plane_wirepitch[Nplanes]/D"); // wire pitch of each plane
    fGeoTree->Branch("plane_wireangle"  , &fPlane_wireangle, "plane_wireangle[Nplanes]/D"); // wire pitch of each plane
    fGeoTree->Branch("plane_wires" , &fPlane_wires, "Plane_wires[Nplanes]/I"); // number of wires in each plane

    fGeoTree->Branch("Nchannels" , &fNchannels);  // number of total channels
    fGeoTree->Branch("NOpDet" , &fNOpDet);  // number of total optical detectors (PMTs)


    // init Event TTree
    TDirectory* subDir2 = fOutFile->mkdir("Event");
    subDir2->cd();
    fEventTree = new TTree("Sim", "Event Tree from Simulation");
    fEventTree->Branch("eventNo", &fEvent);
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);

    fEventTree->Branch("raw_Nhit", &fRaw_Nhit);  // number of hit channels above threshold
    fEventTree->Branch("raw_channelId" , &fRaw_channelId, "raw_channelId[raw_Nhit]/I"); // hit channel id; size == raw_Nhit
    fEventTree->Branch("raw_baseline" , &fRaw_baseline, "raw_baseline[raw_Nhit]/I"); // hit channel baseline; size == raw_Nhit
    fEventTree->Branch("raw_charge" , &fRaw_charge, "raw_charge[raw_Nhit]/I"); // hit channel id; size == raw_Nhit
    fEventTree->Branch("raw_time" , &fRaw_time, "raw_time[raw_Nhit]/I"); // hit channel id; size == raw_Nhit
    fEventTree->Branch("raw_wfADC", &fRaw_wfADC);  // raw waveform adc of each channel
    fEventTree->Branch("raw_wfTDC", &fRaw_wfTDC);  // raw waveform tdc of each channel

    fEventTree->Branch("calib_Nhit", &fCalib_Nhit);  // number of hit channels above threshold
    fEventTree->Branch("calib_channelId" , &fCalib_channelId, "calib_channelId[calib_Nhit]/I"); // hit channel id; size == calib_Nhit
    fEventTree->Branch("calib_wfADC", &fCalib_wfADC);  // calib waveform adc of each channel
    fEventTree->Branch("calib_wfTDC", &fCalib_wfTDC);  // calib waveform tdc of each channel


    fEventTree->Branch("mc_Ntrack", &fMC_Ntrack);  // number of tracks in MC
    fEventTree->Branch("mc_id", &fMC_id, "mc_id[mc_Ntrack]/I");  // track id; size == mc_Ntrack
    fEventTree->Branch("mc_pdg", &fMC_pdg, "mc_id[mc_Ntrack]/I");  // track particle pdg; size == mc_Ntrack
    fEventTree->Branch("mc_process", &fMC_process, "mc_process[mc_Ntrack]/I");  // track generation process code; size == mc_Ntrack    
    fEventTree->Branch("mc_mother", &fMC_mother, "mc_mother[mc_Ntrack]/I");  // mother id of this track; size == mc_Ntrack
    fEventTree->Branch("mc_daughters", &fMC_daughters);  // daughters id of this track; vector
    fEventTree->Branch("mc_startXYZT", &fMC_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
    fEventTree->Branch("mc_endXYZT", &fMC_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
    fEventTree->Branch("mc_startMomentum", &fMC_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
    fEventTree->Branch("mc_endMomentum", &fMC_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack

    fEventTree->Branch("mc_isnu", &mc_isnu);  
    fEventTree->Branch("mc_nGeniePrimaries", &mc_nGeniePrimaries);  
    fEventTree->Branch("mc_nu_pdg", &mc_nu_pdg);
    fEventTree->Branch("mc_nu_ccnc", &mc_nu_ccnc);
    fEventTree->Branch("mc_nu_mode", &mc_nu_mode);
    fEventTree->Branch("mc_nu_intType", &mc_nu_intType);
    fEventTree->Branch("mc_nu_target", &mc_nu_target);
    fEventTree->Branch("mc_hitnuc", &mc_hitnuc);
    fEventTree->Branch("mc_hitquark", &mc_hitquark);
    fEventTree->Branch("mc_nu_Q2", &mc_nu_Q2);
    fEventTree->Branch("mc_nu_W", &mc_nu_W);
    fEventTree->Branch("mc_nu_X", &mc_nu_X);
    fEventTree->Branch("mc_nu_Y", &mc_nu_Y);
    fEventTree->Branch("mc_nu_Pt", &mc_nu_Pt);
    fEventTree->Branch("mc_nu_Theta", &mc_nu_Theta);
    fEventTree->Branch("mc_nu_pos", &mc_nu_pos, "mc_nu_pos[4]/F");
    fEventTree->Branch("mc_nu_mom", &mc_nu_mom, "mc_nu_mom[4]/F");

    fEventTree->Branch("no_hits", &no_hits);  //number of hits
    fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[no_hits]/I");  // channel ID
    fEventTree->Branch("hit_plane", &hit_plane, "hit_plane[no_hits]/I");  // channel plane
    fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[no_hits]/F");  // peak time
    fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[no_hits]/F");  // charge (area)

    for (size_t i=0; i<fTrackModuleLabel.size(); i++) {
        TString name = fTrackModuleLabel[i];
        fEventTree->Branch((name+"_nTrack"    ).Data(), &trk_nTrack[i]);
        fEventTree->Branch((name+"_nHit"      ).Data(), &trk_nHit[i]);
        fEventTree->Branch((name+"_length"    ).Data(), &trk_length[i]);
        fEventTree->Branch((name+"_start_xyz" ).Data(), &trk_start_xyz[i]);
        fEventTree->Branch((name+"_end_xyz"   ).Data(), &trk_end_xyz[i]);
        fEventTree->Branch((name+"_start_dxyz").Data(), &trk_start_dxyz[i]);
        fEventTree->Branch((name+"_end_dxyz"  ).Data(), &trk_end_dxyz[i]);

        fEventTree->Branch((name+"_points_x" ).Data(), &trk_points_x[i]);
        fEventTree->Branch((name+"_points_y" ).Data(), &trk_points_y[i]);
        fEventTree->Branch((name+"_points_z" ).Data(), &trk_points_z[i]);

        for (int j=0; j<3; j++) {
            TString caloName;
            caloName.Form("%s_calo%i_KE", name.Data(), j);
            fEventTree->Branch(caloName.Data(), &trk_calo_KE[i][j]);
            caloName.Form("%s_calo%i_range", name.Data(), j);
            fEventTree->Branch(caloName.Data(), &trk_calo_range[i][j]);
            caloName.Form("%s_calo%i_nHit", name.Data(), j);
            fEventTree->Branch(caloName.Data(), &trk_calo_nHit[i][j]);
            caloName.Form("%s_calo%i_dedx", name.Data(), j);
            fEventTree->Branch(caloName.Data(), &trk_calo_dedx[i][j]);
            caloName.Form("%s_calo%i_dqdx", name.Data(), j);
            fEventTree->Branch(caloName.Data(), &trk_calo_dqdx[i][j]);
            caloName.Form("%s_calo%i_resRange", name.Data(), j);
            fEventTree->Branch(caloName.Data(), &trk_calo_resRange[i][j]);
        }
        
    }

    gDirectory = tmpDir;

}

//-----------------------------------------------------------------------
void CTree::beginJob()
{
 
    fTPC_x = fGeom->DetHalfWidth(0)*2;
    fTPC_y = fGeom->DetHalfHeight(0)*2;
    fTPC_z = fGeom->DetLength(0);

    fNplanes = fGeom->Nplanes();
    for (int i=0; i<fNplanes; i++) {
        fPlane_type[i] = fGeom->Plane(i).SignalType();
        fPlane_view[i] = fGeom->Plane(i).View();
        fPlane_wirepitch[i] = fGeom->WirePitch(fPlane_view[i]);
        fPlane_wireangle[i] = fGeom->WireAngleToVertical(fGeom->Plane(i).View());
        fPlane_wires[i] = fGeom->Nwires(i);
    }
    fNchannels = fGeom->Nchannels();
    fNOpDet = fGeom->NOpDet();
    printGeometry();

    // Write fGeoTree to Disk (once)
    fGeoTree->Fill();
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Detector");
    fGeoTree->Write();
    gDirectory = tmpDir;
    

    // Save Channel Map to text file.
    if (fSaveChannelWireMap) {
        saveChannelWireMap();
        saveWireGeometry();
        saveChannelWireGeometry();
    }

}

//-----------------------------------------------------------------------
void CTree::saveChannelWireMap()
{
    ofstream mapfile;
    mapfile.open("ChannelWireMap.txt");
    mapfile << "# total channels: " << fNchannels << "\n";
    mapfile << "channel\tplane\twire\n"; 
    for (int i=0; i<fNchannels; i++) {
        geo::WireID wid = fGeom->ChannelToWire(i).at(0);
        mapfile << i << "\t" << wid.Plane << "\t" << wid.Wire << "\n";
    }
    mapfile.close();
}

//-----------------------------------------------------------------------
void CTree::saveWireGeometry()
{
    ofstream wirefile;
    wirefile.open("WireGeometry.txt");
    int cstat = 0;
    int tpc = 0;
    double xyzStart[3];
    double xyzEnd[3];

    wirefile << "plane\twire\tsx\tsy\tsz\tex\tey\tez\n"; 
    for (int plane=0; plane<3; plane++) {
        int Nwires = fGeom->Nwires(plane, tpc);
        for (int wire=0; wire<Nwires; wire++) {
            fGeom->WireEndPoints(cstat, tpc, plane, wire, xyzStart, xyzEnd);
            wirefile << plane << "\t" << wire << "\t";
            for (int i=0; i<3; i++) {
                wirefile << xyzStart[i] << "\t";
            }
            for (int i=0; i<3; i++) {
                wirefile << xyzEnd[i] << "\t";
            }
            wirefile << "\n";
        }
    }
    wirefile.close();
}

//-----------------------------------------------------------------------
void CTree::saveChannelWireGeometry()
{
    ofstream out;
    out.open("ChannelWireGeometry.txt");
    int cstat = 0;
    int tpc = 0;
    double xyzStart[3];
    double xyzEnd[3];
    out << "channel\tplane\twire\tsx\tsy\tsz\tex\tey\tez\n"; 
    for (int i=0; i<fNchannels; i++) {
        geo::WireID wid = fGeom->ChannelToWire(i).at(0);
        int plane = wid.Plane;
        int wire = wid.Wire;

        fGeom->WireEndPoints(cstat, tpc, plane, wire, xyzStart, xyzEnd);

        out << i << "\t" << plane << "\t" << wire << "\t";
        for (int i=0; i<3; i++) {
            out << xyzStart[i] << "\t";
        }
        for (int i=0; i<3; i++) {
            out << xyzEnd[i] << "\t";
        }
        out << "\n";
    }
    out.close();
}

//-----------------------------------------------------------------------
void CTree::endJob()
{
    // Write fEventTree to fEventTree
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Event");
    fEventTree->Write();
    gDirectory = tmpDir;

    fOutFile->Close();
}

//-----------------------------------------------------------------------
void CTree::beginRun(const art::Run& /*run*/)
{
  mf::LogInfo("CTree") << "begin run";
}

//-----------------------------------------------------------------------
void CTree::analyze( const art::Event& event )
{
    reset();
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    processMC(event);
    processRaw(event);
    processCalib(event);
    processHits(event);
    processTrack(fTrackModuleLabel[0], event, 0);
    // printEvent();
    fEventTree->Fill();
}

//-----------------------------------------------------------------------
void CTree::reset()
{
    fRaw_wfADC.clear();
    fRaw_wfTDC.clear();
    fCalib_wfADC.clear();
    fCalib_wfTDC.clear();
    for (int i=0; i<MAX_CHANNEL; i++) {
        fRaw_channelId[i] = 0;
        fRaw_baseline[i] = 0;
        fRaw_charge[i] = 0;
        fRaw_time[i] = 0;
        fCalib_channelId[i] = 0;
    }

    for (int i=0; i<MAX_TRACKS; i++) {
        fMC_id[i] = 0;
        fMC_pdg[i] = 0;
        fMC_mother[i] = 0;
        for (int j=0; j<4; j++) {
            fMC_startXYZT[i][j]      = 0;
            fMC_endXYZT[i][j]        = 0;
            fMC_startMomentum[i][j] = 0;
            fMC_endMomentum[i][j]   = 0;
        }
    }
    fMC_daughters.clear();

    for (int i=0; i<MAX_HITS; i++) {
        hit_channel[i] = 0;
        hit_plane[i] = 0;
        hit_peakT[i] = 0;
        hit_charge[i] = 0;
    }

    mc_isnu = 0;
    mc_nGeniePrimaries = -1;
    mc_nu_pdg = -1;
    mc_nu_ccnc = -1;
    mc_nu_mode = -1;
    mc_nu_intType = -1;
    mc_nu_target = -1;
    mc_hitnuc = -1;
    mc_hitquark = -1;
    mc_nu_Q2 = -1;
    mc_nu_W = -1;
    mc_nu_X = -1;
    mc_nu_Y = -1;
    mc_nu_Pt = -1;
    mc_nu_Theta = -1;
    for (int i=0; i<4; i++) {
        mc_nu_pos[i] = 0;
        mc_nu_mom[i] = 0;
    }

    for (size_t i=0; i<MAX_TRACKER; i++) {
        trk_nTrack[i] = -1;
        trk_nHit[i].clear();
        trk_length[i].clear();
        trk_start_xyz[i].clear();
        trk_end_xyz[i].clear();
        trk_start_dxyz[i].clear();
        trk_end_dxyz[i].clear();

        trk_points_x[i].clear();
        trk_points_y[i].clear();
        trk_points_z[i].clear();
        for (size_t j=0; j<3; j++) {
            trk_calo_KE[i][j].clear();
            trk_calo_range[i][j].clear();
            trk_calo_nHit[i][j].clear();
            trk_calo_dedx[i][j].clear();
            trk_calo_dqdx[i][j].clear();
            trk_calo_resRange[i][j].clear();
        }
    }

}

//-----------------------------------------------------------------------
void CTree::processRaw( const art::Event& event )
{
    art::Handle< std::vector<raw::RawDigit> > rawdigit;
    if (! event.getByLabel(fRawDigitLabel, rawdigit)) return;
    std::vector< art::Ptr<raw::RawDigit> >  rawhits;
    art::fill_ptr_vector(rawhits, rawdigit);

    //loop through all RawDigits (over entire channels)
    fRaw_Nhit = 0;
    for (auto const& hit : rawhits) {      
        int chanId = hit->Channel();
        int nSamples = hit->Samples();
        // int pedstal = hit->GetPedestal(); // should be 0 as of 5/13

        std::vector<short> uncompressed(nSamples);
        raw::Uncompress(hit->ADCs(), uncompressed, hit->Compression());


        // determine the baseline of each wire
        int baseline = 0;
        int nBaselineSample = 100;  // use first 100 samples for baseline
        for (int i=0; i<nBaselineSample; i++) { baseline += uncompressed[i]; }
        baseline /= nBaselineSample;
        // cout << chanId << "baseline: " << baseline << endl;

        // uncompressed size is 3200*3 samples per waveform
        short thresh = 2; // threshold set to 2 adc;
        bool isHit = false;
        for (auto const& adc : uncompressed) {
            if (adc > baseline + thresh) {
                isHit = true;
                break;
            }
        }
        if (!isHit) continue; // skip empty channels

        int id = fRaw_Nhit;
        fRaw_channelId[id] = chanId;
        fRaw_baseline[id] = baseline;

        vector<int> wfADC;
        vector<int> wfTDC;
        int nSavedSamples = 0;
        bool hasTime = false;
        for (int i=0; i<nSamples; i++) {
            short adc = uncompressed[i] - baseline;
            if (fabs(adc) > thresh) {
                // cout << i << "," << adc << " | ";
                nSavedSamples++;
                wfADC.push_back(int(adc));
                wfTDC.push_back(i);
                fRaw_charge[id] += adc;
                if (!hasTime) {
                    fRaw_time[id] = i;
                    hasTime = true;
                }
            }
        }

        fRaw_wfADC.push_back(wfADC);
        fRaw_wfTDC.push_back(wfTDC);
        fRaw_Nhit++;
        // cout 
        // << "\n channelID: " << fRaw_channelId[id] 
        // << "\n charge: " << fRaw_charge[id] 
        // << "\n time: " << fRaw_time[id] 
        // << "\n nSamples: " << nSamples
        // << "\n pedestal: " << pedstal
        // << "\n nSavedSamples: " << nSavedSamples
        // << endl;

    }
}

//-----------------------------------------------------------------------
void CTree::processCalib( const art::Event& event )
{

    art::Handle< std::vector<recob::Wire> > wires_handle;
    if (! event.getByLabel(fCalibLabel, wires_handle)) return;
    std::vector< art::Ptr<recob::Wire> >  wires;
    art::fill_ptr_vector(wires, wires_handle);

    // wires size should == Nchannels == 1992; (no hit channel has a flat 0-waveform)
    // cout << "\n wires size: " << wires.size() << endl;

    fCalib_Nhit = 0;
    for (auto const& wire : wires) {      
        std::vector<float> calibwf = wire->Signal(); 
        int chanId = wire->Channel();
        int nSamples = calibwf.size();
        int pedstal = 0;

        short thresh = pedstal + 2; // threshold set to 2 adc;
        bool isHit = false;
        for (auto const& adc : calibwf) {
            if (adc > thresh) {
                isHit = true;
                break;
            }
        }
        if (!isHit) continue; // skip empty channels

        int id = fCalib_Nhit;
        fCalib_channelId[id] = chanId;

        vector<int> wfADC;
        vector<int> wfTDC;
        int nSavedSamples = 0;
        for (int i=0; i<nSamples; i++) {
            int adc = int(calibwf[i]);
            if (fabs(adc) > thresh) {
                // cout << i << "," << adc << " | ";
                nSavedSamples++;
                wfADC.push_back(adc);
                wfTDC.push_back(i);
            }
        }
        // cout << endl;

        fCalib_wfADC.push_back(wfADC);
        fCalib_wfTDC.push_back(wfTDC);
        fCalib_Nhit++;
    }

}

//-----------------------------------------------------------------------
void CTree::processHits( const art::Event& event )
{
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    if (! event.getByLabel(fHitsModuleLabel, hitListHandle)) return;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    art::fill_ptr_vector(hitlist, hitListHandle);

    no_hits = hitlist.size();
    if (no_hits>MAX_HITS) {
        cout << "WARNING:: Event has " << no_hits
             << " hits, MAX_HITS = " << MAX_HITS;
    }

    for (int i=0; i<no_hits; i++) {
        art::Ptr<recob::Hit> hit = hitlist[i];
        hit_channel[i] = hit->Channel();
        hit_plane[i] = hit->WireID().Plane;
        hit_charge[i] = hit->Integral();
        hit_peakT[i] = hit->PeakTime(); 
        if (i==MAX_HITS) break;       
    }
}


//-----------------------------------------------------------------------
void CTree::processTrack(const string& trackLabel,  const art::Event& event, int trackerIndex)
{
    // cout << "track: " << trackLabel << endl;

    art::Handle< std::vector<recob::Track> >  trackListHandle;
    if (! event.getByLabel(trackLabel, trackListHandle)) return;
    std::vector< art::Ptr<recob::Track> > tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);

    trk_nTrack[trackerIndex] = tracklist.size();
    // cout << "nTrack: " << trk_nTrack[trackIndex] << endl;

    art::FindMany<anab::ParticleID> fmpid(trackListHandle, event, trackLabel+"pid");
    art::FindMany<anab::Calorimetry> fmcal(trackListHandle, event, trackLabel+"calo");
    
    for (int i=0; i<trk_nTrack[trackerIndex]; i++) {
        art::Ptr<recob::Track> track = tracklist[i];

        int nHit = track->NumberTrajectoryPoints();
        trk_nHit[trackerIndex].push_back(nHit);
        trk_length[trackerIndex].push_back(TrackLength(*track));

        // cout << "TrackID: " << track->ID() << endl;
        // cout << "NumberTrajectoryPoints: " << track->NumberTrajectoryPoints() << endl;


        vector<double> start_xyz, end_xyz, start_dxyz, end_dxyz;
        start_xyz.push_back(track->Vertex().X());
        start_xyz.push_back(track->Vertex().Y());
        start_xyz.push_back(track->Vertex().Z());
        end_xyz.push_back(track->End().X());
        end_xyz.push_back(track->End().Y());
        end_xyz.push_back(track->End().Z());
        start_dxyz.push_back(track->VertexDirection().X());
        start_dxyz.push_back(track->VertexDirection().Y());
        start_dxyz.push_back(track->VertexDirection().Z());
        end_dxyz.push_back(track->EndDirection().X());
        end_dxyz.push_back(track->EndDirection().Y());
        end_dxyz.push_back(track->EndDirection().Z());
        trk_start_xyz[trackerIndex].push_back(start_xyz);
        trk_end_xyz[trackerIndex].push_back(end_xyz);
        trk_start_dxyz[trackerIndex].push_back(start_dxyz);
        trk_end_dxyz[trackerIndex].push_back(end_dxyz);

        vector<double> points_x, points_y, points_z;
        for (int iPt=0; iPt<nHit; iPt++) {
            TVector3 point = track->LocationAtPoint(iPt);
            points_x.push_back(point.X());
            points_y.push_back(point.Y());
            points_z.push_back(point.Z());
        }
        trk_points_x[trackerIndex].push_back(points_x);
        trk_points_y[trackerIndex].push_back(points_y);
        trk_points_z[trackerIndex].push_back(points_z);

        // cout << "NumberFitMomentum: " << track->NumberFitMomentum() << endl;
        // cout << "start pos: " << track->Vertex().X() << endl;
        // cout << "start dir: " << track->VertexDirection().X() << endl;
        // cout << "start mom: " << track->VertexMomentum() << endl;

        // cout << "end pos: " << track->End().X() << endl;
        // cout << "end dir: " << track->EndDirection().X() << endl;
        // cout << "end mom: " << track->EndMomentum() << endl;

        cout << "track length: " << TrackLength(*track) << endl;

        // add pid info associated with the track
        if (fmpid.isValid()) {
            std::vector<const anab::ParticleID*> pids = fmpid.at(i);
            cout << "pid size: " << pids.size() << endl;
            const anab::ParticleID* pid = pids.at(0);
            cout << "pid pdg: " << pid->Pdg() << endl;
            cout << "pid minChi2: " << pid->MinChi2() << endl;
            cout << "pid Chi2Proton: " << pid->Chi2Proton() << endl;
            cout << "pid Chi2Kaon: " << pid->Chi2Kaon() << endl;
            cout << "pid Chi2Pion: " << pid->Chi2Pion() << endl;
            cout << "pid Chi2Muon: " << pid->Chi2Muon() << endl;
            cout << "pid PIDA: " << pid->PIDA() << endl;
        }
        else {
            cout << "pid: not found" << endl;
        }

        // add calo info associated with the track
        if (fmcal.isValid()){
            std::vector<const anab::Calorimetry*> calos = fmcal.at(i);

            if (calos.size() != 3) {
                cout << "ERROR: calos size: " << calos.size() << ", should be 3!" << endl;
            }

            for (int iCalo=0; iCalo<3; iCalo++) {
                const anab::Calorimetry* calo = calos[i];
                trk_calo_KE[trackerIndex][iCalo].push_back(calo->KineticEnergy());
                trk_calo_range[trackerIndex][iCalo].push_back(calo->Range());
                trk_calo_nHit[trackerIndex][iCalo].push_back(calo->dEdx().size());

                trk_calo_dedx[trackerIndex][iCalo].push_back(calo->dEdx());
                trk_calo_dqdx[trackerIndex][iCalo].push_back(calo->dQdx());
                trk_calo_resRange[trackerIndex][iCalo].push_back(calo->ResidualRange());

                // cout << "calo " << i << endl;
                // cout << "ke: " <<  calo->KineticEnergy() << endl;
                // cout << "range: " << calo->Range() << endl;
                // cout << "pitchc: " << calo->TrkPitchC() << endl;
                // cout << "hit size: " << calo->dEdx().size() << endl;
            }
        }
        else {
            cout <<"calo: not found" << endl;
        }


        cout << endl;

    }

}

//-----------------------------------------------------------------------
void CTree::processMC( const art::Event& event )
{
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    if (! event.getByLabel("largeant", particleHandle)) return;
    std::vector< art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);

    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel("largeant", simChannelHandle);    

    fMC_Ntrack = particles.size();
    if (fMC_Ntrack > MAX_TRACKS) {
        cout << "WARNING:: # tracks " << fMC_Ntrack << " exceeded MAX_TRACKS " << MAX_TRACKS << endl;
    }
    int i=0; // track index
    for (auto const& particle: particles ) {
        fMC_process[i] = processMap[particle->Process()];
        if (fMC_process[i] == 0) cout << "unknown process: " << particle->Process() << endl;
        fMC_id[i] = particle->TrackId();
        fMC_pdg[i] = particle->PdgCode();
        fMC_mother[i] = particle->Mother();
        int Ndaughters = particle->NumberDaughters();
        vector<int> daughters;
        for (int i=0; i<Ndaughters; i++) {
            daughters.push_back(particle->Daughter(i));
        }
        fMC_daughters.push_back(daughters);
        size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
        int last = numberTrajectoryPoints - 1;
        const TLorentzVector& positionStart = particle->Position(0);
        const TLorentzVector& positionEnd   = particle->Position(last);
        const TLorentzVector& momentumStart = particle->Momentum(0);
        const TLorentzVector& momentumEnd   = particle->Momentum(last);
        positionStart.GetXYZT(fMC_startXYZT[i]);
        positionEnd.GetXYZT(fMC_endXYZT[i]);
        momentumStart.GetXYZT(fMC_startMomentum[i]);
        momentumEnd.GetXYZT(fMC_endMomentum[i]);
        i++;
        if (i==MAX_TRACKS) break;
    } // particle loop done 

    // Generator Info 
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    event.getByLabel("generator",mctruthListHandle);
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mctruthListHandle);
    art::Ptr<simb::MCTruth> mctruth;

    if (mclist.size()>0) {
        // cout << "MCTruth size: " << mclist.size() << endl;
        mctruth = mclist.at(0);
        if (mctruth->NeutrinoSet()) {
            simb::MCNeutrino nu = mctruth->GetNeutrino();
            mc_isnu = 1;
            mc_nGeniePrimaries = mctruth->NParticles();
            mc_nu_pdg = nu.Nu().PdgCode();
            mc_nu_ccnc = nu.CCNC();
            mc_nu_mode = nu.Mode();
            mc_nu_intType = nu.InteractionType();
            mc_nu_target = nu.Target();
            mc_hitnuc = nu.HitNuc();
            mc_hitquark = nu.HitQuark();
            mc_nu_Q2 = nu.QSqr();
            mc_nu_W = nu.W();
            mc_nu_X = nu.X();
            mc_nu_Y = nu.Y();
            mc_nu_Pt = nu.Pt();
            mc_nu_Theta = nu.Theta();

            const TLorentzVector& position = nu.Nu().Position(0);
            const TLorentzVector& momentum = nu.Nu().Momentum(0);
            position.GetXYZT(mc_nu_pos);
            momentum.GetXYZT(mc_nu_mom);
        }
    }
    
}

//-----------------------------------------------------------------------
void CTree::printGeometry()
{
    cout << endl;
    cout << "====== microboone geometry ======" << endl;
    cout << "# cryostats: " << fGeom->Ncryostats() << endl;
    cout << "# TPC: " << fGeom->NTPC() << endl;  
    cout << "\tTPC xyz: " << fTPC_x << ", " << fTPC_y << ", " << fTPC_z << endl;
    cout << "# planes: " << fNplanes << endl;
    for (int i=0; i<fNplanes; i++) {
        cout 
            << "\tplane " << i 
            << "( type: " << fPlane_type[i]
            << ", view: " << fPlane_view[i]
            << ", wirepitch: " << fPlane_wirepitch[i]
            << ", wire angle: " << fPlane_wireangle[i]*180/TMath::Pi()
            << ", wires: " << fPlane_wires[i]
            << ")" << endl;
    }
    cout << "# channels: " << fNchannels << endl;
    cout << "# OpDet: " << fGeom->NOpDet() << endl;
    cout << "# AuxDetectors: " << fGeom->NAuxDets() << endl;
    cout << "====== microboone geometry ======" << endl;
    cout << endl;
}

//-----------------------------------------------------------------------
void CTree::printEvent()
{
    cout << " Event: " << fEvent << endl;
    cout << "   Run: " << fRun << endl;
    cout << "SubRun: " << fSubRun << endl;

    cout << "Raw Hit Channels: " << fRaw_Nhit << endl;
    cout << "      Ntracks:" << fMC_Ntrack;
    for (int i=0; i<fMC_Ntrack; i++) {
        cout << "\n              id: " << fMC_id[i];
        cout << "\n             pdg: " << fMC_pdg[i];
        cout << "\n          mother: " << fMC_mother[i];
        cout << "\n      Ndaughters: " << fMC_daughters.at(i).size();
        cout << "\n      start XYZT: (" << fMC_startXYZT[i][0] << ", " << fMC_startXYZT[i][1] << ", " << fMC_startXYZT[i][2] << ", " << fMC_startXYZT[i][3] << ")";
        cout << "\n        end XYZT: (" << fMC_endXYZT[i][0] << ", " << fMC_endXYZT[i][1] << ", " << fMC_endXYZT[i][2] << ", " << fMC_endXYZT[i][3] << ")";
        cout << "\n  start momentum: (" << fMC_startMomentum[i][0] << ", " << fMC_startMomentum[i][1] << ", " << fMC_startMomentum[i][2] << ", " << fMC_startMomentum[i][3] << ")";
        cout << "\n    end momentum: (" << fMC_endMomentum[i][0] << ", " << fMC_endMomentum[i][1] << ", " << fMC_endMomentum[i][2] << ", " << fMC_endMomentum[i][3] << ")";

        cout << endl;
    }

    // cout << "Number of Hits Found: " << no_hits << endl;
    // // for (int i=0; i<no_hits; i++) {
    // //     cout << hit_channel[i] << ", ";
    // //     cout << hit_charge[i] << ", ";
    // //     cout << hit_peakT[i] << ", ";
    // //     cout << endl;
    // // }

}

//-----------------------------------------------------------------------
void CTree::InitProcessMap()
{
    processMap["unknown"]              = 0;
    processMap["primary"]              = 1;
    processMap["compt"]                = 2;
    processMap["phot"]                 = 3;
    processMap["annihil"]              = 4;
    processMap["eIoni"]                = 5;
    processMap["eBrem"]                = 6;
    processMap["conv"]                 = 7;
    processMap["muIoni"]               = 8;
    processMap["muMinusCaptureAtRest"] = 9;
    processMap["NeutronInelastic"]     = 10;
    processMap["nCapture"]             = 11;
    processMap["hadElastic"]           = 12;
}

//-----------------------------------------------------------------------
double CTree::TrackLength(const recob::Track& track)
{
    double result = 0.;   
    TVector3 disp = track.LocationAtPoint(0);   
    int n = track.NumberTrajectoryPoints();      
    for(int i = 1; i < n; ++i) {    
        const TVector3& pos = track.LocationAtPoint(i);    
        //double momentum = track.MomentumAtPoint(i);    
        //std::cout<<"\n"<<i<<"\t"<<momentum<<"\n";    
        disp -= pos;    
        result += disp.Mag();    
        disp = pos;
    }    
    return result;
}


//-----------------------------------------------------------------------
DEFINE_ART_MODULE(CTree)
} // namespace microboone


#endif
