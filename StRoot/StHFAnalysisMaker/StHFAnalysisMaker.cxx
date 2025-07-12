#include "StHFAnalysisMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "CutsConfig.h" 
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include <cmath>
#include <iostream>

ClassImp(StHFAnalysisMaker)

StHFAnalysisMaker::StHFAnalysisMaker(const char* name):StMaker(name) {}
StHFAnalysisMaker::~StHFAnalysisMaker(){}

Int_t StHFAnalysisMaker::Init(){
    if(!mPicoDstMaker){LOG_ERROR<<"HFAnalysis: set PicoDstMaker first"<<endm;return kStFatal;}
    // book hists
    hJPsiMass = new TH1F("hJPsiMass","e^{+}e^{-} inv mass;M [GeV]",120,2.0,4.0);
    hJPsiPtY  = new TH2F("hJPsiPtY","J/#psi pT vs y;pT;y",100,0,10,60,-HFCuts::PID::nSigmaE,HFCuts::PID::nSigmaE);
    hD0Mass   = new TH1F("hD0Mass","K#pi inv mass;M [GeV]",120,1.6,2.1);
    hD0PtY    = new TH2F("hD0PtY","D^{0} pT vs y;pT;y",100,0,10,60,-HFCuts::PID::nSigmaE,HFCuts::PID::nSigmaE);
    hNPEPt    = new TH1F("hNPEPt","NPE pT;pT",100,0,10);
    hEoverPvsP= new TH2F("hEoverPvsP","E/p vs p;p;E/p",100,0,10,100,0,2);

    // additional observables
    hPhiVsEP_JPsi = new TH2F("hPhiVsEP_JPsi","J/#psi #phi-#Psi_{2} vs p_{T};p_{T};#phi-#Psi_{2}",100,0,10,72,-TMath::Pi(),TMath::Pi());
    hPhiVsEP_D0   = new TH2F("hPhiVsEP_D0","D^{0} #phi-#Psi_{2} vs p_{T};p_{T};#phi-#Psi_{2}",100,0,10,72,-TMath::Pi(),TMath::Pi());
    hED0_DeltaPhi = new TH1F("hED0_DeltaPhi","e-D^{0} #Delta#phi;#Delta#phi",72,-TMath::Pi(),TMath::Pi());
    hEOPInclusive = new TH1F("hEOPInclusive","Inclusive e E/p;E/p",100,0,2);
    hEffMap_JPsi  = new TH2F("hEffMap_JPsi","J/#psi counts (proxy for eff);p_{T};y",100,0,10,60,-3,3);
    hEffMap_D0    = new TH2F("hEffMap_D0","D^{0} counts (proxy for eff);p_{T};y",100,0,10,60,-3,3);
    hRefMultVz    = new TH2F("hRefMultVz","gRefMult vs Vz;Vz (cm);gRefMult",120,-60,60,100,0,1000);
    return kStOK;
}

bool StHFAnalysisMaker::passEventCuts(){
    return HFCuts::passEventCuts(mPicoDstMaker->picoDst()->event());
}

bool StHFAnalysisMaker::goodTrack(const StPicoTrack* t){
    if(!t) return false;
    if(t->nHitsFit()<HFCuts::Track::nHitsFitMin) return false;
    if(t->nHitsFit()/(float)t->nHitsPoss()<HFCuts::Track::nHitsFracMin) return false;
    auto vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
    if(t->gDCA(TVector3(vtx.x(),vtx.y(),vtx.z())).Mag() > HFCuts::Track::dcaMax) return false;
    return true;
}







Int_t StHFAnalysisMaker::Make(){
    if(!passEventCuts()) return kStOK;

    // Single loop for particle identification and kinematic calculation
    auto pico = mPicoDstMaker->picoDst();
    const int nTr = pico->numberOfTracks();
    std::vector<Particle> electrons, kaons, pions;

    for(int i=0; i<nTr; ++i){
        auto t = pico->track(i);
        if(!goodTrack(t)) continue;

        TVector3 pMom = t->pMom();

        if(std::fabs(t->nSigmaElectron()) < HFCuts::PID::nSigmaE && pMom.Perp() > HFCuts::Track::ptMin){
            electrons.push_back({pMom, std::sqrt(pMom.Mag2() + 0.000511*0.000511), t->charge()});
        }
        if(std::fabs(t->nSigmaKaon()) < HFCuts::PID::nSigmaK){
            kaons.push_back({pMom, std::sqrt(pMom.Mag2() + 0.493677*0.493677), t->charge()});
        }
        if(std::fabs(t->nSigmaPion()) < HFCuts::PID::nSigmaPi){
            pions.push_back({pMom, std::sqrt(pMom.Mag2() + 0.13957*0.13957), t->charge()});
        }
    }

    // event-level QA
    auto evt = mPicoDstMaker->picoDst()->event();
    if(evt) hRefMultVz->Fill(evt->primaryVertex().z(),evt->grefMult());

    // J/Psi analysis
    for(size_t ia=0; ia<electrons.size(); ++ia){
        const auto& e1 = electrons[ia];
        for(size_t ib=ia+1; ib<electrons.size(); ++ib){
            const auto& e2 = electrons[ib];
            if(e1.charge * e2.charge >=0) continue;
            TVector3 p = e1.pMom + e2.pMom;
            double e = e1.energy + e2.energy;
            double m=std::sqrt(e*e-p.Mag2()); double pt=p.Perp(); double y=0.5*std::log((e+p.Z())/(e-p.Z()+1e-6));
            hJPsiMass->Fill(m); hJPsiPtY->Fill(pt,y); hEffMap_JPsi->Fill(pt,y);
            double dphi = std::atan2(p.Y(),p.X()); // proxy phi
            hPhiVsEP_JPsi->Fill(pt,dphi);
        }
    }

    // D0 analysis
    for(const auto& kaon : kaons){
        for(const auto& pion : pions){
            if(kaon.charge * pion.charge >= 0) continue;
            TVector3 q = kaon.pMom + pion.pMom;
            double e = kaon.energy + pion.energy;
            double m=std::sqrt(e*e-q.Mag2()); double pt=q.Perp(); double y=0.5*std::log((e+q.Z())/(e-q.Z()+1e-6));
            hD0Mass->Fill(m);
            hD0PtY->Fill(pt,y);
            hEffMap_D0->Fill(pt,y);

            // Mass window cut for e-D0 correlation
            if(m > HFCuts::Mass::D0_M_min && m < HFCuts::Mass::D0_M_max){
                for(const auto& el : electrons){
                    double dphi = std::fabs(el.pMom.Phi() - q.Phi());
                    if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                    hED0_DeltaPhi->Fill(dphi);
                }
            }
        }
    }

    // HFE analysis
    for(const auto& el : electrons){
        double p = el.pMom.Mag();
        // TODO: retrieve BEMC pid traits for E/p. Skip if traits unavailable.
        hEoverPvsP->Fill(p,-1);
    }
    return kStOK;
}
    
Int_t StHFAnalysisMaker::Finish(){
    if(mOutFile.empty()){LOG_ERROR<<"Output file name is empty"<<endm; return kStFatal;}
    TFile *f = new TFile(mOutFile.c_str(),"RECREATE","HFAnalysisOutput",9);
    if(!f||f->IsZombie()){LOG_ERROR<<"Cannot create output file "<<mOutFile<<endm; return kStFatal;}
    TH1* hists[] = {hJPsiMass,hD0Mass,hNPEPt,hEOPInclusive};
    for(auto h: hists) if(h) h->Write();
    TH2* h2s[] = {hJPsiPtY,hD0PtY,hEoverPvsP,hPhiVsEP_JPsi,hPhiVsEP_D0,hEffMap_JPsi,hEffMap_D0,hRefMultVz};
    for(auto h:h2s) if(h) h->Write();
    f->Close();
    return kStOK;
}
