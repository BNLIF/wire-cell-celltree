// Read data from MC raw files and convert it into ROOT tree
// Chao Zhang (chao@bnl.gov) 5/13/2014

#ifndef CELLTREE_MODULE
#define CELLTREE_MODULE

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
#include "TClonesArray.h"
#include "TH1F.h"

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

class CellTree : public art::EDAnalyzer {
public:

    explicit CellTree(fhicl::ParameterSet const& pset);
    virtual ~CellTree();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void initOutput();
    void printEvent();


    void processCalib(const art::Event& evt);
    void processSimChannel(const art::Event& evt);
    void reset();


private:

    // the parameters we'll read from the .fcl
    std::string fCalibLabel;
    std::string fOutFileName;
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

    // art::ServiceHandle<geo::Geometry> fGeom;
    // // art::ServiceHandle<util::LArProperties> larp;

    TFile *fOutFile;
    TTree *fEventTree;

    // Event Tree Leafs
    int fEvent;
    int fRun;
    int fSubRun;


    int fCalib_nChannel;
    // int fCalib_channelId[MAX_CHANNEL];  // hit channel id; size == fCalib_Nhit
    // // FIXEME:: cannot save e.g std::vector<std::vector<float> > in ttree
    std::vector<int> fCalib_channelId;
    // std::vector<std::vector<float> > fCalib_wf;
    TClonesArray *fCalib_wf;
    // std::vector<std::vector<int> > fCalib_wfTDC;

    int fSIMIDE_size;
    vector<int> fSIMIDE_channelIdY;
    vector<int> fSIMIDE_trackId;
    vector<unsigned short> fSIMIDE_tdc;
    vector<float> fSIMIDE_x;
    vector<float> fSIMIDE_y;
    vector<float> fSIMIDE_z;
    vector<float> fSIMIDE_numElectrons;

}; // class CellTree


//-----------------------------------------------------------------------
CellTree::CellTree(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
    initOutput();
}

//-----------------------------------------------------------------------
CellTree::~CellTree()
{
}

//-----------------------------------------------------------------------
void CellTree::reconfigure(fhicl::ParameterSet const& p){
    fCalibLabel = p.get<std::string>("CalibLabel");
    fOutFileName = p.get<std::string>("outFile");
}

//-----------------------------------------------------------------------
void CellTree::initOutput()
{
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");

    TNamed version("version", "2.0");
    version.Write();

    // init Event TTree
    TDirectory* subDir = fOutFile->mkdir("Event");
    subDir->cd();
    fEventTree = new TTree("Sim", "Event Tree from Simulation");
    fEventTree->Branch("eventNo", &fEvent);
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);

    fEventTree->Branch("calib_nChannel", &fCalib_nChannel);  // number of hit channels above threshold
    fEventTree->Branch("calib_channelId" , &fCalib_channelId); // hit channel id; size == calib_Nhit

    fCalib_wf = new TClonesArray("TH1F");
    fEventTree->Branch("calib_wf", &fCalib_wf, 256000, 0);  // calib waveform adc of each channel
    // fCalib_wf->BypassStreamer();
    // fEventTree->Branch("calib_wfTDC", &fCalib_wfTDC);  // calib waveform tdc of each channel

    fEventTree->Branch("simide_size", &fSIMIDE_size);  // size of stored sim:IDE
    fEventTree->Branch("simide_channelIdY", &fSIMIDE_channelIdY);
    fEventTree->Branch("simide_trackId", &fSIMIDE_trackId);
    fEventTree->Branch("simide_tdc", &fSIMIDE_tdc);
    fEventTree->Branch("simide_x", &fSIMIDE_x);
    fEventTree->Branch("simide_y", &fSIMIDE_y);
    fEventTree->Branch("simide_z", &fSIMIDE_z);
    fEventTree->Branch("simide_numElectrons", &fSIMIDE_numElectrons);

    gDirectory = tmpDir;

}

//-----------------------------------------------------------------------
void CellTree::beginJob()
{


}


//-----------------------------------------------------------------------
void CellTree::endJob()
{
    // Write fEventTree to fEventTree
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Event");

    fEventTree->Write();

    gDirectory = tmpDir;

    fOutFile->Close();
}

//-----------------------------------------------------------------------
void CellTree::beginRun(const art::Run& /*run*/)
{
  mf::LogInfo("CellTree") << "begin run";
}

//-----------------------------------------------------------------------
void CellTree::analyze( const art::Event& event )
{
    reset();
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    processCalib(event);
    processSimChannel(event);

    printEvent();
    fEventTree->Fill();
}

//-----------------------------------------------------------------------
void CellTree::reset()
{

    fCalib_channelId.clear();
    fCalib_wf->Clear();

    fSIMIDE_channelIdY.clear();
    fSIMIDE_trackId.clear();
    fSIMIDE_tdc.clear();
    fSIMIDE_x.clear();
    fSIMIDE_y.clear();
    fSIMIDE_z.clear();
    fSIMIDE_numElectrons.clear();

}

//-----------------------------------------------------------------------
void CellTree::processCalib( const art::Event& event )
{

    art::Handle< std::vector<recob::Wire> > wires_handle;
    if (! event.getByLabel(fCalibLabel, wires_handle)) {
        cout << "WARNING: no label " << fCalibLabel << endl;
        return;
    }
    std::vector< art::Ptr<recob::Wire> >  wires;
    art::fill_ptr_vector(wires, wires_handle);

    // wires size should == Nchannels == 1992; (no hit channel has a flat 0-waveform)
    // cout << "\n wires size: " << wires.size() << endl;
    fCalib_nChannel = wires.size();

    int i=0;
    for (auto const& wire: wires) {
        std::vector<float> calibwf = wire->Signal();
        int chanId = wire->Channel();
        fCalib_channelId.push_back(chanId);
        TH1F *h = new((*fCalib_wf)[i]) TH1F("","", 9600, 0, 9600);
        for (int j=1; j<=9600; j++) {
            h->SetBinContent(j, calibwf[j]);
        }
        // fCalib_wf.push_back(calibwf);
        // cout << chanId << ", " << nSamples << endl;
        i++;
    }

}

//-----------------------------------------------------------------------
void CellTree::processSimChannel( const art::Event& event )
{
    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel("largeant", simChannelHandle);
    cout << "total simChannel: " << (*simChannelHandle).size() << endl;
    fSIMIDE_size = 0;
    for ( auto const& channel : (*simChannelHandle) ) {
        auto channelNumber = channel.Channel();
        // cout << channelNumber << endl;
        // if (! (fGeometry->SignalType( channelNumber ) == geo::kCollection) ) {
        //     continue;
        // }
        auto const& timeSlices = channel.TDCIDEMap();
        for ( auto const& timeSlice : timeSlices ) {
            auto const& energyDeposits = timeSlice.second;
            for ( auto const& energyDeposit : energyDeposits ) {
                fSIMIDE_size++;
                fSIMIDE_channelIdY.push_back(channelNumber);
                fSIMIDE_tdc.push_back(timeSlice.first);
                fSIMIDE_trackId.push_back(energyDeposit.trackID);
                fSIMIDE_x.push_back(energyDeposit.x);
                fSIMIDE_y.push_back(energyDeposit.y);
                fSIMIDE_z.push_back(energyDeposit.z);
                fSIMIDE_numElectrons.push_back(energyDeposit.numElectrons);
                // cout << channelNumber << ": " << energyDeposit.trackID << ": " << timeSlice.first << ": "
                //      << energyDeposit.x << ", " << energyDeposit.y << ", " << energyDeposit.z << ", "
                //      << energyDeposit.numElectrons << endl;
            }
        }

    }
    cout << "total IDEs: " << fSIMIDE_size << endl;


}


//-----------------------------------------------------------------------
void CellTree::printEvent()
{
    cout << " Event: " << fEvent << endl;
    cout << "   Run: " << fRun << endl;
    cout << "SubRun: " << fSubRun << endl;

}


//-----------------------------------------------------------------------
DEFINE_ART_MODULE(CellTree)
} // namespace microboone


#endif
