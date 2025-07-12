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
    // dielectron like/unlike-sign spectra
    hMee_LSneg  = new TH1F("hMee_LSneg","e^{-}e^{-} mass;M [GeV]",120,0,4);
    hMee_LSpos  = new TH1F("hMee_LSpos","e^{+}e^{+} mass;M [GeV]",120,0,4);
    hMee_ULS    = new TH1F("hMee_ULS","e^{+}e^{-} mass;M [GeV]",120,0,4);
    hMeePt_LSneg= new TH2F("hMeePt_LSneg","e^{-}e^{-} mass vs p_{T};M;p_{T}",120,0,4,100,0,10);
    hMeePt_LSpos= new TH2F("hMeePt_LSpos","e^{+}e^{+} mass vs p_{T};M;p_{T}",120,0,4,100,0,10);
    hMeePt_ULS  = new TH2F("hMeePt_ULS","e^{+}e^{-} mass vs p_{T};M;p_{T}",120,0,4,100,0,10);
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
    const size_t nE = mElectrons.size();
    for(size_t ia=0; ia<nE; ++ia){
        const auto* e1 = mElectrons[ia];
        if(e1->pMom().Perp()<HFCuts::Track::ptMin) continue;
        for(size_t ib=ia+1; ib<nE; ++ib){
            const auto* e2 = mElectrons[ib];
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
    const std::vector<const StPicoTrack*>& kp = mKplus;
    const std::vector<const StPicoTrack*>& kpbar = mKminus;
    const std::vector<const StPicoTrack*>& pi = mPiplus;
    const std::vector<const StPicoTrack*>& piBar = mPiminus;
    const auto& electrons = mElectrons;
    auto comb=[&](const std::vector<const StPicoTrack*>& K,const std::vector<const StPicoTrack*>& P){
        for(size_t ik=0; ik<K.size(); ++ik){
            const auto* k = K[ik];
            for(size_t jp=0; jp<P.size(); ++jp){
                const auto* p = P[jp];
                if(k==p) continue;
                TVector3 pk=k->pMom(), pp=p->pMom(), q=pk+pp;
            double e=std::sqrt(pk.Mag2()+0.493677*0.493677)+std::sqrt(pp.Mag2()+0.13957*0.13957);
            double m=std::sqrt(e*e-q.Mag2()); double pt=q.Perp(); double y=0.5*std::log((e+q.Z())/(e-q.Z()+1e-6));
            hD0Mass->Fill(m); hD0PtY->Fill(pt,y); hEffMap_D0->Fill(pt,y);
                // compute e-D0 delta phi with electrons from earlier
                for(const auto* et : electrons){
                    double dphi = std::fabs(et->pMom().Phi() - q.Phi());
                    if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                    hED0_DeltaPhi->Fill(dphi);
                }
            } 
        } 
    };
    comb(kpbar,pi); comb(kp,piBar);
}

void StHFAnalysisMaker::runHFE(){
    for(const auto* t : mElectrons){
        double p = t->pMom().Mag();
        hEoverPvsP->Fill(p,-1);
    }
}

Int_t StHFAnalysisMaker::Make(){
    if(!passEventCuts()) return kStOK;
    // build cached track lists
    mElectrons.clear(); mKplus.clear(); mKminus.clear(); mPiplus.clear(); mPiminus.clear();
    auto picoEvt = mPicoDstMaker->picoDst();
    const int nTrTot = picoEvt->numberOfTracks();
    for(int i=0;i<nTrTot;++i){
        auto t = picoEvt->track(i);
        if(!goodTrack(t)) continue;
        const float nSigE = std::fabs(t->nSigmaElectron());
        const float nSigK = std::fabs(t->nSigmaKaon());
        const float nSigPi= std::fabs(t->nSigmaPion());
        if(nSigE < HFCuts::PID::nSigmaE) mElectrons.push_back(t);
        if(nSigK < HFCuts::PID::nSigmaK) ((t->charge()>0)? mKplus : mKminus).push_back(t);
        if(nSigPi< HFCuts::PID::nSigmaPi) ((t->charge()>0)? mPiplus: mPiminus).push_back(t);
    }
    // event-level QA
    auto evt = picoEvt->event();
    if(evt) hRefMultVz->Fill(evt->primaryVertex().z(),evt->grefMult());

    runJPsi(); runD0(); runHFE(); runDielectronPairs();
    return kStOK;
}

Int_t StHFAnalysisMaker::Finish(){
    if(mOutFile.empty()){LOG_ERROR<<"Output file name is empty"<<endm; return kStFatal;}
    TFile *f = new TFile(mOutFile.c_str(),"RECREATE","HFAnalysisOutput",9);
    if(!f||f->IsZombie()){LOG_ERROR<<"Cannot create output file "<<mOutFile<<endm; return kStFatal;}
    TH1* hists[] = {hJPsiMass,hD0Mass,hNPEPt,hEOPInclusive,hMee_LSneg,hMee_LSpos,hMee_ULS};
    for(auto h: hists) if(h) h->Write();
    TH2* h2s[] = {hJPsiPtY,hD0PtY,hEoverPvsP,hPhiVsEP_JPsi,hPhiVsEP_D0,hEffMap_JPsi,hEffMap_D0,hRefMultVz,
                       hMeePt_LSneg,hMeePt_LSpos,hMeePt_ULS};
    for(auto h:h2s) if(h) h->Write();
    f->Close();
    return kStOK;
}

// --- Dielectron pairing ---
void StHFAnalysisMaker::runDielectronPairs(){
    struct ETrack{const StPicoTrack* tr; TLorentzVector lv;};
    std::vector<ETrack> plus, minus; plus.reserve(mElectrons.size()); minus.reserve(mElectrons.size());
    const double me=0.000511; // GeV
    for(const auto* t: mElectrons){
        TVector3 p = t->pMom();
        TLorentzVector lv(p, std::sqrt(p.Mag2()+me*me));
        (t->charge()>0 ? plus : minus).push_back({t, lv});
    }
    auto pairLoop=[&](const std::vector<ETrack>& A,const std::vector<ETrack>& B,
                      TH1* hM,TH2* hMPt,bool sameList){
        for(size_t i=0;i<A.size();++i){
            const TLorentzVector& l1 = A[i].lv;
            size_t jStart = sameList ? i+1 : 0;
            for(size_t j=jStart;j<B.size();++j){
                const TLorentzVector& l2 = B[j].lv;
                TLorentzVector pr = l1 + l2;
                if(hM)   hM->Fill(pr.M());
                if(hMPt) hMPt->Fill(pr.M(), pr.Pt());
            }
        }
    };
    pairLoop(minus, minus, hMee_LSneg , hMeePt_LSneg, true);
    pairLoop(plus , plus , hMee_LSpos , hMeePt_LSpos, true);
    pairLoop(plus , minus, hMee_ULS   , hMeePt_ULS  , false);
}

