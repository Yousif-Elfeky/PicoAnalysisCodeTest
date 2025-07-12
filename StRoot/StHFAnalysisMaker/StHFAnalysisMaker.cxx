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
    hRefMultVz = new TH2F("hRefMultVz","gRefMult vs Vz;Vz (cm);gRefMult",120,-60,60,100,0,1000);

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

void StHFAnalysisMaker::runJPsi(){
    auto pico = mPicoDstMaker->picoDst();
    std::vector<int> els;
    const int nTr = pico->numberOfTracks();
    for(int i=0;i<nTr;++i){
        auto tr=pico->track(i);
        if(!goodTrack(tr)) continue;
        if(std::fabs(tr->nSigmaElectron())>HFCuts::PID::nSigmaE) continue;
        if(tr->pMom().Perp()<HFCuts::Track::ptMin) continue;
        els.push_back(i);
    }
    for(size_t ia=0; ia<els.size(); ++ia){
        auto e1=pico->track(els[ia]);
        for(size_t ib=ia+1; ib<els.size(); ++ib){
            auto e2=pico->track(els[ib]);
            if(e1->charge()*e2->charge()>=0) continue;
            TVector3 p1=e1->pMom(), p2=e2->pMom();
            TVector3 p=p1+p2; double e=std::sqrt(p1.Mag2()+0.000511*0.000511)+std::sqrt(p2.Mag2()+0.000511*0.000511);
            double m=std::sqrt(e*e-p.Mag2()); double pt=p.Perp(); double y=0.5*std::log((e+p.Z())/(e-p.Z()+1e-6));
            hJPsiMass->Fill(m); hJPsiPtY->Fill(pt,y); hEffMap_JPsi->Fill(pt,y);
            double dphi = std::atan2(p.Y(),p.X()); // proxy phi
            hPhiVsEP_JPsi->Fill(pt,dphi);
        }
    }
}

void StHFAnalysisMaker::runD0(){
    auto pico=mPicoDstMaker->picoDst();
    std::vector<const StPicoTrack*> kp,kpbar,pi,piBar;
    const int nTr=pico->numberOfTracks();
    for(int i=0;i<nTr;++i){
        auto t=pico->track(i);
        if(!goodTrack(t)) continue;
        if(std::fabs(t->nSigmaKaon())<HFCuts::PID::nSigmaK) (t->charge()>0?kp:kpbar).push_back(t);
        if(std::fabs(t->nSigmaPion())<HFCuts::PID::nSigmaPi) (t->charge()>0?pi:piBar).push_back(t);
    }
    auto comb=[&](const std::vector<const StPicoTrack*>& K,const std::vector<const StPicoTrack*>& P){
        for(auto k:K){for(auto p:P){ if(k==p) continue; TVector3 pk=k->pMom(), pp=p->pMom(), q=pk+pp;
            double e=std::sqrt(pk.Mag2()+0.493677*0.493677)+std::sqrt(pp.Mag2()+0.13957*0.13957);
            double m=std::sqrt(e*e-q.Mag2()); double pt=q.Perp(); double y=0.5*std::log((e+q.Z())/(e-q.Z()+1e-6));
            hD0Mass->Fill(m); hD0PtY->Fill(pt,y); hEffMap_D0->Fill(pt,y);
                // compute e-D0 delta phi with electrons from earlier
                for(int ie=0; ie<pico->numberOfTracks(); ++ie){ 
                    auto et=pico->track(ie); 
                    if(!goodTrack(et)) continue; 
                    if(std::fabs(et->nSigmaElectron())>HFCuts::PID::nSigmaE) continue;
                    double dphi=std::fabs(et->pMom().Phi()-q.Phi()); 
                    if(dphi> TMath::Pi()) dphi=2*TMath::Pi()-dphi; 
                    hED0_DeltaPhi->Fill(dphi);
                } 
            } 
        } 
    };
    comb(kpbar,pi); comb(kp,piBar);
}

void StHFAnalysisMaker::runHFE(){
    auto pico=mPicoDstMaker->picoDst();
    const int nTr=pico->numberOfTracks();
    for(int i=0;i<nTr;++i){auto t=pico->track(i); if(!goodTrack(t)) continue; if(std::fabs(t->nSigmaElectron())>HFCuts::PID::nSigmaE) continue;
        double p=t->pMom().Mag();
        // TODO: retrieve BEMC pid traits for E/p. Skip if traits unavailable.
        hEoverPvsP->Fill(p,-1); }
}

Int_t StHFAnalysisMaker::Make(){
    if(!passEventCuts()) return kStOK;
    // event-level QA
    auto evt = mPicoDstMaker->picoDst()->event();
    if(evt) hRefMultVz->Fill(evt->primaryVertex().z(),evt->grefMult());

    runJPsi(); runD0(); runHFE();
    return kStOK;
}

Int_t StHFAnalysisMaker::Finish(){
    TFile* f=TFile::Open(mOutFile.c_str(),"RECREATE"); if(!f||f->IsZombie()){LOG_ERROR<<"Cannot write output"<<endm; return kStFatal;}
    TH1* hists[] = {hJPsiMass,hD0Mass,hNPEPt,hEOPInclusive};
    for(auto h: hists) if(h) h->Write();
    TH2* h2s[] = {hJPsiPtY,hD0PtY,hEoverPvsP,hPhiVsEP_JPsi,hPhiVsEP_D0,hEffMap_JPsi,hEffMap_D0,hRefMultVz};
    for(auto h:h2s) if(h) h->Write();
    f->Close();
    return kStOK;
}
