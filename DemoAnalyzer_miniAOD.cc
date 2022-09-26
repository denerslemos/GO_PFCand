// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dener De Souza Lemos
//         Created:  Wed, 12 Feb 2020 13:05:50 GMT
//
//



// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//v0 candidates
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

//
// class declaration
//
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include "TVector3.h"
#include <vector>
#include <map>
#include <functional>

#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include <TLorentzVector.h>

#define pimass 0.1396


//
// CMSSW user include files
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include "DataFormats/HeavyIonEvent/interface/HFFilterInfo.h" //this line is needed to access the HF Filters

//#include "trackingEfficiency2018PbPb.h"



using reco::TrackCollection;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );
      
      void initHistos(const edm::Service<TFileService> & fs);  

      Double_t GetQ(const TLorentzVector &p1, const TLorentzVector &p2);

      const TLorentzVector InvertXYVector( TLorentzVector &vec);

      bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2, TH2F* cosdpt);

      void hbt(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<int> GoodTrackCharge, Double_t cc, THnSparse* h_ss, THnSparse* h_os, THnSparse* h_ss_rot, THnSparse* h_os_rot, TH2F* cosdpt);
      
      void MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<int>> ev_GoodTrackCharge, THnSparse* h_ss, THnSparse* h_os); 

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
//      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracksToken_;
      edm::EDGetTokenT<edm::ValueMap<float>> chi2Map_;
      edm::EDGetTokenT<reco::HFFilterInfo> HFfilters_;
      edm::EDGetTokenT<reco::Centrality> CentralityTag_;
      edm::EDGetTokenT<int> CentralityBinTag_; 
  	  
      double vtx_x, vtx_y, vtx_z, vtx_rho, vtx_xError, vtx_yError, vtx_zError;  //primary vertex
  
      int cent;
	  
      TH1F* Cent;

      TH1F* Vtx_z;
      TH1F* Vtx_x;
      TH1F* Vtx_y;
      TH1F* Vtx_z_err;
      TH1F* Vtx_x_err;
      TH1F* Vtx_y_err;

      TH2F* pT_reco;
      TH2F* pT_reco_corr;
      TH2F* pT_reco_corr_pos;
      TH2F* pT_reco_corr_neg;

      TH2F* eta_reco;
      TH2F* eta_reco_corr;
      TH2F* eta_reco_corr_neg;
      TH2F* eta_reco_corr_pos;

      TH2F* phi_reco;
      TH2F* phi_reco_corr;
      TH2F* phi_reco_corr_neg;
      TH2F* phi_reco_corr_pos;

      TH1F* NHits_reco;
      TH1F* Algo_reco;
      TH1F* dzSig_reco;
      TH1F* dxySig_reco;
      TH1F* pTres_reco;
      TH1F* chi2_reco;

      TH1F* NHits_reco_corr;
      TH1F* Algo_reco_corr;
      TH1F* dzSig_reco_corr;
      TH1F* dxySig_reco_corr;
      TH1F* pTres_reco_corr;
      TH1F* chi2_reco_corr; 

      THnSparseF *hist_SS;
      THnSparseF *hist_OS;
      THnSparseF *hist_SS_mix;
      THnSparseF *hist_OS_mix;

      THnSparseF *hist_SS_corr;
      THnSparseF *hist_OS_corr;
      THnSparseF *hist_SS_mix_corr;
      THnSparseF *hist_OS_mix_corr;

      TH2F* cosineDpt;
      TH2F* cosineDpt_corr;

      std::vector<double> ev_z_vtx;
      std::vector<int>    ev_cc_vec;
      std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec;
      std::vector<std::vector<int>>  ev_GoodTrackCharge_vec;

      std::vector<double> ev_z_vtx_corr;
      std::vector<int>    ev_cc_vec_corr;
      std::vector<std::vector<TLorentzVector>>  ev_GoodTrackFourVector_vec_corr;
      std::vector<std::vector<int>>  ev_GoodTrackCharge_vec_corr;

      int numMinHFTower4;

      int ev=0;
      int evpass=0,evfail=0;

      int evpassx=0,evfailx=0;
      int evpassy=0,evfaily=0;

      int ccc = 0;

      int evNumber;

      TH1F* tower_phi;
      TH1F* tower_eta;
      TH1F* tower_pt;
      TH1F* tower_et;
      TH1F* tower_energy;
      TH1F* tower_energy_lead;
      TH1F* tower_ethad;
      TH1F* tower_etem;
      TH1F* tower_etsum;

      TH1F* Ntowers;
      TH1F* Ntower_TH2GeV_plus;
      TH1F* Ntower_TH4GeV_plus;
      TH1F* Ntower_TH2GeV_minus;
      TH1F* Ntower_TH4GeV_minus;
      TH1F* MinNtower_TH2GeV;
      TH1F* MinNtower_TH4GeV;

/*

TrkEff2018PbPb trkEf_ =  TrkEff2018PbPb("general", false, "");

TrkEff2018PbPb trkEff_plus =  TrkEff2018PbPb("generalMB+", false, "");

TrkEff2018PbPb trkEff_minus =  TrkEff2018PbPb("generalMB-", false, "");

*/

};

DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
 :
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  tracksToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks"))),
  chi2Map_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("chi2src"))),
  HFfilters_(consumes<reco::HFFilterInfo>(iConfig.getParameter<edm::InputTag>("HFfilters"))),
  CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"))),
  CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc")))
{
   //now do what ever initialization is needed
}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

//Vertex selection

  evNumber      = -9999;

  int evn = 0;
  int evnx = 0;

  edm::Handle<std::vector<reco::Vertex> > vertexCollection;
  iEvent.getByToken(vertexToken_,vertexCollection);
  std::vector<reco::Vertex> vtx_sorted = *vertexCollection;
  std::sort( vtx_sorted.begin(), vtx_sorted.end(), DemoAnalyzer::vtxSort );
  if(vtx_sorted.size() == 0) return;
  vtx_x = (double)vtx_sorted.begin()->position().x();// VertexX_->Fill(vtx_x);
  vtx_y = (double)vtx_sorted.begin()->position().y();// VertexY_->Fill(vtx_y);
  vtx_z = (double)vtx_sorted.begin()->position().z(); 
  vtx_rho= (double)vtx_sorted.begin()->position().Rho(); 
  vtx_xError = (double)vtx_sorted.begin()->xError();
  vtx_yError = (double)vtx_sorted.begin()->yError();
  vtx_zError = (double)vtx_sorted.begin()->zError();
  math::XYZPoint vtxx(vtx_x,vtx_y,vtx_z);

//  if(fabs(vtx_z) > 15.)return;

  Vtx_z->Fill(vtx_z);
  Vtx_x->Fill(vtx_x);
  Vtx_y->Fill(vtx_y);
  Vtx_z_err->Fill(vtx_zError);
  Vtx_x_err->Fill(vtx_xError);
  Vtx_y_err->Fill(vtx_yError);

  edm::Handle<int> cbin_;
  iEvent.getByToken(CentralityBinTag_,cbin_);
  cent = *cbin_;
  Cent->Fill(cent);
    
//  Handle<TrackCollection> tracks;
  edm::Handle<edm::View<pat::PackedCandidate>> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<edm::ValueMap<float>> chi2Map;
  iEvent.getByToken(chi2Map_,chi2Map);

  edm::Handle<reco::HFFilterInfo> HFfilter;
  iEvent.getByToken(HFfilters_, HFfilter);
  numMinHFTower4 = HFfilter->numMinHFTowers4;

  evNumber     = iEvent.id().event();

  ev=ev+1;
  if(HFfilter->numMinHFTowers4 < 2){evfail=evfail+1; cout << "Event fail: " <<  evNumber << endl; evn = evNumber;}else{evpass=evpass+1;}
  if(numMinHFTower4 < 2){evfailx=evfailx+1;}else{evpassx=evpassx+1;}

//  if(numMinHFTower4 < 2)return;

  std::vector<TLorentzVector> GoodTrackFourVector;
  std::vector<int> GoodTrackCharge;
  std::vector<TLorentzVector> GoodTrackFourVector_corr;
  std::vector<int> GoodTrackCharge_corr;


//  for(TrackCollection::const_iterator iter_tk = tracks->begin(); iter_tk != tracks->end();++iter_tk) {
  unsigned short int count = 0;
  unsigned short int count_TH2_plus = 0;
  unsigned short int count_TH2_minus = 0;
  unsigned short int count_TH4_plus = 0;
  unsigned short int count_TH4_minus = 0;
  float leading_energy = 0;
    //Loop over packedPF candidates
    for(unsigned it = 0; it<tracks->size(); ++it){
      const pat::PackedCandidate & tk = (*tracks)[it];
      
      //pdgId needed: 1 for hadronic and 2 for EM particles in HF
      if(tk.pdgId()==1 || tk.pdgId()==2){
          const bool eta_plus = (tk.eta() > 3.0) && (tk.eta() < 6.0);
          const bool eta_minus = (tk.eta() < -3.0) && (tk.eta() > -6.0);
          if(tk.et()<0.0)continue;
          if(eta_plus || eta_minus){
         	 tower_eta->Fill(tk.eta());
           	 tower_phi->Fill(tk.phi());
	  	 tower_pt->Fill(tk.pt());
          	 tower_et->Fill(tk.et());
	  	 tower_energy->Fill(tk.energy());
		 if(tk.energy()>=leading_energy){leading_energy = tk.energy();}
	  	 if(tk.pdgId()==1)tower_ethad->Fill(tk.et());
          	 if(tk.pdgId()==2)tower_etem->Fill(tk.et());
	  	 count = count+1;
          }
          if(eta_minus){
		count_TH2_minus += tk.energy() >= 2.0 ? 1 : 0;
                count_TH4_minus += tk.energy() >= 4.0 ? 1 : 0;
	  } else if(eta_plus){
		count_TH2_plus += tk.energy() >= 2.0 ? 1 : 0;
                count_TH4_plus += tk.energy() >= 4.0 ? 1 : 0;
	  }
      }

    }

    tower_energy_lead->Fill(leading_energy);

    Ntowers->Fill(count);
    Ntower_TH2GeV_plus->Fill(count_TH2_plus);
    Ntower_TH4GeV_plus->Fill(count_TH4_plus);
    Ntower_TH2GeV_minus->Fill(count_TH2_minus);
    Ntower_TH4GeV_minus->Fill(count_TH4_minus);

    MinNtower_TH2GeV->Fill(std::min(count_TH2_plus, count_TH2_minus));
    MinNtower_TH4GeV->Fill(std::min(count_TH4_plus, count_TH4_minus));

    if(std::min(count_TH4_plus, count_TH4_minus) < 2){evfaily=evfaily+1; cout << "Event fail: " <<  evNumber << endl; evnx = evNumber;}else{evpassy=evpassy+1;}

if(evn==evnx && evn!=0) ccc = ccc + 1;

/*
      if(!tk.hasTrackDetails()) continue;

      reco::Track const& iter_tk = tk.pseudoTrack();

      double aux_tk_dz_vtx = (double)iter_tk.dz(vtxx);
      double aux_tk_dzError_vtx  = (double)sqrt(iter_tk.dzError()*iter_tk.dzError()+vtx_zError*vtx_zError);
      double aux_tk_dxy_vtx = (double)iter_tk.dxy(vtxx);
      double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk.dxyError()*iter_tk.dxyError()+vtx_xError*vtx_yError);
      double chi2ndof = (double) (*chi2Map)[tracks->ptrAt(it)];
      if(iter_tk.pt()<=0.5)continue;
//      if(!iter_tk->quality(reco::TrackBase::highPurity))continue; // not needed in PF Candidates
      if(fabs(iter_tk.eta())>2.4)continue;
      const reco::HitPattern& hit_pattern = iter_tk.hitPattern();
      //make all plots and vectors without tracking selection

      TLorentzVector pvector;
      pvector.SetXYZM(iter_tk.px(),iter_tk.py(),iter_tk.pz(),pimass);
//      GoodTrackFourVector.push_back(pvector);
//      GoodTrackCharge.push_back(iter_tk->charge());

      pT_reco->Fill(iter_tk.pt(),cent);
      eta_reco->Fill(iter_tk.eta(),cent);
      phi_reco->Fill(iter_tk.phi(),cent);

      NHits_reco->Fill(iter_tk.numberOfValidHits());
      Algo_reco->Fill(iter_tk.algo());
      dzSig_reco->Fill(aux_tk_dz_vtx/aux_tk_dzError_vtx);
      dxySig_reco->Fill(aux_tk_dxy_vtx/aux_tk_dxyError_vtx);
      pTres_reco->Fill(fabs(iter_tk.ptError())/iter_tk.pt());
      chi2_reco->Fill(chi2ndof/hit_pattern.trackerLayersWithMeasurement());

      if(fabs(iter_tk.ptError())/iter_tk.pt()>0.1)continue;
      if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
      if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue;
      if(hit_pattern.pixelLayersWithMeasurement()>0){GoodTrackFourVector.push_back(pvector); GoodTrackCharge.push_back(iter_tk.charge());}
      if((chi2ndof/hit_pattern.trackerLayersWithMeasurement())>0.18)continue;
      if(iter_tk.numberOfValidHits()<11)continue;
      if(hit_pattern.pixelLayersWithMeasurement()>0){GoodTrackFourVector_corr.push_back(pvector); GoodTrackCharge_corr.push_back(iter_tk.charge());}

      pT_reco_corr->Fill(iter_tk.pt(),cent);
      eta_reco_corr->Fill(iter_tk.eta(),cent);
      phi_reco_corr->Fill(iter_tk.phi(),cent);

      NHits_reco_corr->Fill(iter_tk.numberOfValidHits());
      Algo_reco_corr->Fill(iter_tk.algo());
      dzSig_reco_corr->Fill(aux_tk_dz_vtx/aux_tk_dzError_vtx);
      dxySig_reco_corr->Fill(aux_tk_dxy_vtx/aux_tk_dxyError_vtx);
      pTres_reco_corr->Fill(fabs(iter_tk.ptError())/iter_tk.pt());
      chi2_reco_corr->Fill(chi2ndof/hit_pattern.trackerLayersWithMeasurement());

      if(iter_tk.charge()>0){
      pT_reco_corr_pos->Fill(iter_tk.pt(),cent);
      eta_reco_corr_pos->Fill(iter_tk.eta(),cent);
      phi_reco_corr_pos->Fill(iter_tk.phi(),cent);
      }else if(iter_tk.charge()<0){ 
      pT_reco_corr_neg->Fill(iter_tk.pt(),cent);
      eta_reco_corr_neg->Fill(iter_tk.eta(),cent);
      phi_reco_corr_neg->Fill(iter_tk.phi(),cent);
      }

    }

if(GoodTrackFourVector.size() > 1){

hbt(GoodTrackFourVector, GoodTrackCharge, (Double_t) cent, hist_SS, hist_OS, hist_SS_mix, hist_OS_mix, cosineDpt);

//ev_z_vtx.push_back(vtx_z);
//ev_cc_vec.push_back(cent);
//ev_GoodTrackFourVector_vec.push_back(GoodTrackFourVector);
//ev_GoodTrackCharge_vec.push_back(GoodTrackCharge);

}

if(GoodTrackFourVector_corr.size() > 1){

hbt(GoodTrackFourVector_corr, GoodTrackCharge_corr, (Double_t) cent, hist_SS_corr, hist_OS_corr, hist_SS_mix_corr, hist_OS_mix_corr, cosineDpt_corr);

//ev_z_vtx_corr.push_back(vtx_z);
//ev_cc_vec_corr.push_back(cent);
//ev_GoodTrackFourVector_vec_corr.push_back(GoodTrackFourVector_corr);
//ev_GoodTrackCharge_vec_corr.push_back(GoodTrackCharge_corr);

}
*/

}


// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
  std::cout<<"  This is called once for each job: Begin Job " << std::endl;
  edm::Service<TFileService> fs;
  initHistos(fs);

  cent = -1;
  vtx_x = -999.;
  vtx_y = -999.;
  vtx_z = -999.;
  vtx_rho = -999.;
  vtx_xError = -999.; 
  vtx_yError = -999.;
  vtx_zError = -999.;
  numMinHFTower4 = -1;

}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
std::cout << std::endl;
//std::cout << "Time for mixing" << std::endl;
std::cout << std::endl;

std::cout << "Total of events: " << ev << std::endl;
std::cout << "Passed (standard): " << evpass << std::endl;
std::cout << "Failed (standard): " << evfail << std::endl;
//std::cout << "PassedX: " << evpassx << std::endl;
//std::cout << "FailedX: " << evfailx << std::endl;
std::cout << "Passed (PF candidates): " << evpassy << std::endl;
std::cout << "Failed (PF candidates): " << evfaily << std::endl;

std::cout << "Failed events matched (std vs PF): " << ccc << std::endl;

//without selection
/*
std::cout << "Mixing without selection" << std::endl;
MixEvents(0,9,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "5%" << std::endl;
MixEvents(10,19,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "10%" << std::endl;
MixEvents(20,29,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "15%" << std::endl;
MixEvents(30,39,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "20%" << std::endl;
MixEvents(40,49,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "25%" << std::endl;
MixEvents(50,59,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "30%" << std::endl;
MixEvents(60,69,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "35%" << std::endl;
MixEvents(70,79,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "40%" << std::endl;
MixEvents(80,89,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "45%" << std::endl;
MixEvents(90,99,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "50%" << std::endl;
MixEvents(100,109,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "55%" << std::endl;
MixEvents(110,119,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "60%" << std::endl;
MixEvents(120,129,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "65%" << std::endl;
MixEvents(130,139,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "70%" << std::endl;
MixEvents(140,149,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "75%" << std::endl;
MixEvents(150,159,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "80%" << std::endl;
MixEvents(160,169,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "85%" << std::endl;
MixEvents(170,179,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "90%" << std::endl;
MixEvents(180,189,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "95%" << std::endl;
MixEvents(190,199,10,ev_cc_vec,ev_z_vtx,2.0,ev_GoodTrackFourVector_vec,ev_GoodTrackCharge_vec,hist_SS_mix,hist_OS_mix);
std::cout << "100%" << std::endl;

//with selection
std::cout << "Mixing with selection" << std::endl;
MixEvents(0,9,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "5%" << std::endl;
MixEvents(10,19,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "10%" << std::endl;
MixEvents(20,29,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "15%" << std::endl;
MixEvents(30,39,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "20%" << std::endl;
MixEvents(40,49,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "25%" << std::endl;
MixEvents(50,59,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "30%" << std::endl;
MixEvents(60,69,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "35%" << std::endl;
MixEvents(70,79,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "40%" << std::endl;
MixEvents(80,89,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "45%" << std::endl;
MixEvents(90,99,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "50%" << std::endl;
MixEvents(100,109,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "55%" << std::endl;
MixEvents(110,119,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "60%" << std::endl;
MixEvents(120,129,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "65%" << std::endl;
MixEvents(130,139,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "70%" << std::endl;
MixEvents(140,149,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "75%" << std::endl;
MixEvents(150,159,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "80%" << std::endl;
MixEvents(160,169,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "85%" << std::endl;
MixEvents(170,179,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "90%" << std::endl;
MixEvents(180,189,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "95%" << std::endl;
MixEvents(190,199,10,ev_cc_vec_corr,ev_z_vtx_corr,2.0,ev_GoodTrackFourVector_vec_corr,ev_GoodTrackCharge_vec_corr,hist_SS_mix_corr,hist_OS_mix_corr);
std::cout << "100%" << std::endl;
*/
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

bool DemoAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b ){
   if( a.tracksSize() != b.tracksSize() )
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   else
      return  a.chi2() < b.chi2() ? true : false ;
}

void DemoAnalyzer::initHistos(const edm::Service<TFileService> & fs){

   TFileDirectory Inf = fs->mkdir( "Info" );

   Cent=Inf.make<TH1F>("Cent","",200,0.0,200.0);

   Vtx_z=Inf.make<TH1F>("Vtx_z","",200,-20.0,20.0);
   Vtx_x=Inf.make<TH1F>("Vtx_x","",1000,-1.0,1.0);
   Vtx_y=Inf.make<TH1F>("Vtx_y","",1000,-1.0,1.0);
   Vtx_z_err=Inf.make<TH1F>("Vtx_z_err","",1000,0.0,2.0);
   Vtx_x_err=Inf.make<TH1F>("Vtx_x_err","",1000,0.0,2.0);
   Vtx_y_err=Inf.make<TH1F>("Vtx_y_err","",1000,0.0,2.0);

   pT_reco=Inf.make<TH2F>("pT_reco","",200,0.0,10.0,200,0.0,200.0);
   pT_reco_corr=Inf.make<TH2F>("pT_reco_corr","",200,0.0,10.0,200,0.0,200.0);
   pT_reco_corr_pos=Inf.make<TH2F>("pT_reco_corr_pos","",200,0.0,10.0,200,0.0,200.0);
   pT_reco_corr_neg=Inf.make<TH2F>("pT_reco_corr_neg","",200,0.0,10.0,200,0.0,200.0);

   eta_reco=Inf.make<TH2F>("eta_reco","",100,-2.5,2.5,200,0.0,200.0);
   eta_reco_corr=Inf.make<TH2F>("eta_reco_corr","",100,-2.5,2.5,200,0.0,200.0);
   eta_reco_corr_neg=Inf.make<TH2F>("eta_reco_corr_neg","",100,-2.5,2.5,200,0.0,200.0);
   eta_reco_corr_pos=Inf.make<TH2F>("eta_reco_corr_pos","",100,-2.5,2.5,200,0.0,200.0);

   phi_reco=Inf.make<TH2F>("phi_reco","",100,-3.2,3.2,200,0.0,200.0);
   phi_reco_corr=Inf.make<TH2F>("phi_reco_corr","",100,-3.2,3.2,200,0.0,200.0);
   phi_reco_corr_neg=Inf.make<TH2F>("phi_reco_corr_neg","",100,-3.2,3.2,200,0.0,200.0);
   phi_reco_corr_pos=Inf.make<TH2F>("phi_reco_corr_pos","",100,-3.2,3.2,200,0.0,200.0);

   NHits_reco=Inf.make<TH1F>("NHits_reco","",60,0.0,60.0);
   Algo_reco=Inf.make<TH1F>("Algo_reco","",30,0.0,30.0);
   dzSig_reco=Inf.make<TH1F>("dzSig_reco","",100,-20.0,20.0);
   dxySig_reco=Inf.make<TH1F>("dxySig_reco","",100,-20.0,20.0);
   pTres_reco=Inf.make<TH1F>("pTres_reco","",100,0.0,1.0);
   chi2_reco=Inf.make<TH1F>("chi2_reco","",100,0.0,2.0);

   NHits_reco_corr=Inf.make<TH1F>("NHits_reco_corr","",60,0.0,60.0);
   Algo_reco_corr=Inf.make<TH1F>("Algo_reco_corr","",30,0.0,30.0);
   dzSig_reco_corr=Inf.make<TH1F>("dzSig_reco_corr","",100,-20.0,20.0);
   dxySig_reco_corr=Inf.make<TH1F>("dxySig_reco_corr","",100,-20.0,20.0);
   pTres_reco_corr=Inf.make<TH1F>("pTres_reco_corr","",100,0.0,1.0);
   chi2_reco_corr=Inf.make<TH1F>("chi2_reco_corr","",100,0.0,2.0);


   Int_t bins3D[3]=   {100,15,10};
   Double_t xmin3D[3]={0.,0.,0.};
   Double_t xmax3D[3]={2.,1.5,200.};

   hist_SS=Inf.make<THnSparseF>("hist_SS","hist_SS",3,bins3D,xmin3D,xmax3D);
   hist_OS=Inf.make<THnSparseF>("hist_OS","hist_OS",3,bins3D,xmin3D,xmax3D);
   hist_SS_mix=Inf.make<THnSparseF>("hist_SS_mix","hist_SS_mix",3,bins3D,xmin3D,xmax3D);
   hist_OS_mix=Inf.make<THnSparseF>("hist_OS_mix","hist_OS_mix",3,bins3D,xmin3D,xmax3D);

   hist_SS_corr=Inf.make<THnSparseF>("hist_SS_corr","hist_SS_corr",3,bins3D,xmin3D,xmax3D);
   hist_OS_corr=Inf.make<THnSparseF>("hist_OS_corr","hist_OS_corr",3,bins3D,xmin3D,xmax3D);
   hist_SS_mix_corr=Inf.make<THnSparseF>("hist_SS_mix_corr","hist_SS_mix_corr",3,bins3D,xmin3D,xmax3D);
   hist_OS_mix_corr=Inf.make<THnSparseF>("hist_OS_mix_corr","hist_OS_mix_corr",3,bins3D,xmin3D,xmax3D);
 
   cosineDpt=Inf.make<TH2F>("cosineDpt_reco","",10000,0.999,1.001,1000,-0.1,0.1);
   cosineDpt_corr=Inf.make<TH2F>("cosineDpt_reco_corr","",10000,0.999,1.001,1000,-0.1,0.1);


   tower_pt=Inf.make<TH1F>("tower_pt","",1000,0.0,100.0);
   tower_eta=Inf.make<TH1F>("tower_eta","",200,-6.0,6.0);
   tower_phi=Inf.make<TH1F>("tower_phi","",100,-3.2,3.2);
   tower_et=Inf.make<TH1F>("tower_et","",1000,0.0,100.0);
   tower_energy=Inf.make<TH1F>("tower_energy","",1000,0.0,100.0);
   tower_energy_lead=Inf.make<TH1F>("tower_energy_lead","",1000,0.0,100.0);
   tower_ethad=Inf.make<TH1F>("tower_ethad","",1000,0.0,100.0);
   tower_etem=Inf.make<TH1F>("tower_etem","",1000,0.0,100.0);
   tower_etsum=Inf.make<TH1F>("tower_etsum","",1000,0.0,100.0);

   Ntowers=Inf.make<TH1F>("Ntowers","",1000,0.0,1000.0);
   Ntower_TH2GeV_plus=Inf.make<TH1F>("Ntower_TH2GeV_plus","",1000,0.0,1000.0);
   Ntower_TH4GeV_plus=Inf.make<TH1F>("Ntower_TH4GeV_plus","",1000,0.0,1000.0);
   Ntower_TH2GeV_minus=Inf.make<TH1F>("Ntower_TH2GeV_minus","",1000,0.0,1000.0);
   Ntower_TH4GeV_minus=Inf.make<TH1F>("Ntower_TH4GeV_minus","",1000,0.0,1000.0);
   MinNtower_TH2GeV=Inf.make<TH1F>("MinNtower_TH2GeV","",1000,0.0,1000.0);
   MinNtower_TH4GeV=Inf.make<TH1F>("MinNtower_TH4GeV","",1000,0.0,1000.0);

}

Double_t DemoAnalyzer::GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
  TLorentzVector Sum4V = p1+p2;
  Double_t q = Sum4V.Mag2() - 2*(p1.M()*p1.M() + p2.M()*p2.M());
  return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

const TLorentzVector DemoAnalyzer::InvertXYVector( TLorentzVector &vec){
  TLorentzVector ovec = vec;
  ovec.SetX(-vec.X());
  ovec.SetY(-vec.Y());
  return ovec;
} //DONE

void DemoAnalyzer::hbt(std::vector<TLorentzVector> GoodTrackFourVector, std::vector<int> GoodTrackCharge, Double_t cc, THnSparse* h_ss, THnSparse* h_os, THnSparse* h_ss_rot, THnSparse* h_os_rot, TH2F* cosdpt){
  for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
	int charge1 = GoodTrackCharge[itk1];
    for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
	int charge2 = GoodTrackCharge[itk2];

        Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
    	TLorentzVector psum2 = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
        Double_t kt=(psum2.Pt())/2.;
        Double_t x3D[3]={q,kt,(Double_t)cc};

	Double_t q_rot = GetQ(GoodTrackFourVector[itk1], InvertXYVector(GoodTrackFourVector[itk2]));        
        TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
        Double_t kt_rot = (psum2_rot.Pt())/2.;
        Double_t x3D_rot[3]={q_rot,kt_rot,(Double_t)cc};

	if(splitcomb(GoodTrackFourVector[itk1],  GoodTrackFourVector[itk2], cosdpt))continue;

	if(charge1*charge2>0){h_ss->Fill(x3D);h_ss_rot->Fill(x3D_rot);}else{h_os->Fill(x3D);h_os_rot->Fill(x3D_rot);}
}}}


void DemoAnalyzer::MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> pv_vtx_z, double vzcut,  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector, std::vector<std::vector<int>> ev_GoodTrackCharge, THnSparse* h_ss, THnSparse* h_os){
	int aux_n_evts = (int)ev_ntrkoff.size();
	for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
		std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec = ev_GoodTrackFourVector[nevt1];
		std::vector<int> ev_GoodTrackCharge_nevt1_vec = ev_GoodTrackCharge[nevt1];
		int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
		if(ev_ntrkoff[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt1])continue;
		int takeAssociated = 0;
		for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
			if(ev_ntrkoff[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff[nevt_assoc])continue;
			if(fabs(pv_vtx_z[nevt1] - pv_vtx_z[nevt_assoc]) > vzcut)continue;
			std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec = ev_GoodTrackFourVector[nevt_assoc];
			std::vector<int> ev_GoodTrackCharge_nevt_assoc_vec = ev_GoodTrackCharge[nevt_assoc];
			int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
			takeAssociated++;
			if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.
			for(int imix=0; imix<nMixmult_nevt1; imix++){
				for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
					Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
					TLorentzVector psum2mix = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
					Double_t ktmix=(psum2mix.Pt())/2.;
					Double_t xmix3D[3]={qmix,ktmix,(Double_t)ev_ntrkoff[nevt1]};
					if(ev_GoodTrackCharge_nevt1_vec[imix]*ev_GoodTrackCharge_nevt_assoc_vec[iimix]>0){h_ss->Fill(xmix3D);}else{h_os->Fill(xmix3D);}
				}
			}
		}
	}
}


bool DemoAnalyzer::splitcomb(TLorentzVector &vec1,TLorentzVector &vec2, TH2F* cosdpt){
   bool issplit=false;
   Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());
   Double_t dpt = vec1.Pt() - vec2.Pt();
   cosdpt->Fill(cosa,dpt);
//   if ( (cosa > 0.99999) && (deltapt < 0.02)) { issplit = true;} //same used by Sunil in pPb
   return issplit;
} //DONE


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
edm::ParameterSetDescription desc;
desc.setUnknown();
descriptions.addDefault(desc);
}       


//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
