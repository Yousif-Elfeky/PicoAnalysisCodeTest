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
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TRandom.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "TLorentzVector.h"
#include <cmath>
#include <iostream>
#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
ClassImp(StHFAnalysisMaker)

StHFAnalysisMaker::StHFAnalysisMaker(const char* name):StMaker(name) {}
StHFAnalysisMaker::~StHFAnalysisMaker(){}

Int_t StHFAnalysisMaker::Init(){
    if(!mPicoDstMaker){LOG_ERROR<<"HFAnalysis: set PicoDstMaker first"<<endm;return kStFatal;}
    // book hists
    hJPsiMass = new TH1F("hJPsiMass","e^{+}e^{-} inv mass;M [GeV]",1000,2.0,4.0);
    hJPsiPtY  = new TH2F("hJPsiPtY","J/#psi pT vs y;pT;y",1000,0,10,1000,-HFCuts::PID::nSigmaE,HFCuts::PID::nSigmaE);
    hD0Mass   = new TH1F("hD0Mass","K#pi inv mass;M [GeV]",1000,1.6,2.1);
    hD0PtY    = new TH2F("hD0PtY","D^{0} pT vs y;pT;y",1000,0,10,1000,-HFCuts::PID::nSigmaE,HFCuts::PID::nSigmaE);
    // background mass spectra
    hJPsiBkgMass1 = new TH1F("hJPsiBkgMass1","like-sign e^{+}e^{+} mass;M [GeV]",1000,2.0,4.0);
    hJPsiBkgMass2 = new TH1F("hJPsiBkgMass2","like-sign e^{-}e^{-} mass;M [GeV]",1000,2.0,4.0);
    hD0BkgMass1   = new TH1F("hD0BkgMass1","same-charge K+#pi^+ mass;M [GeV]",1000,1.6,2.1);
    hD0BkgMass2   = new TH1F("hD0BkgMass2","same-charge K-#pi^- mass;M [GeV]",1000,1.6,2.1);

    hNPEPt    = new TH1F("hNPEPt","NPE pT;pT",1000,0,10);
    hEoverPvsP= new TH2F("hEoverPvsP","E/p vs p;p;E/p",1000,0,10,1000,0,2);

    // additional observables
    hPhiVsEP_JPsi = new TH2F("hPhiVsEP_JPsi","J/#psi #phi-#Psi_{2} vs p_{T};p_{T};#phi-#Psi_{2}",1000,0,10,1000,-TMath::Pi(),TMath::Pi());
    hPhiVsEP_D0   = new TH2F("hPhiVsEP_D0","D^{0} #phi-#Psi_{2} vs p_{T};p_{T};#phi-#Psi_{2}",1000,0,10,1000,-TMath::Pi(),TMath::Pi());
    hED0_DeltaPhi = new TH1F("hED0_DeltaPhi","e-D^{0} #Delta#phi;#Delta#phi",1000,-TMath::Pi(),TMath::Pi());
    hEOPInclusive = new TH1F("hEOPInclusive","Inclusive e E/p;E/p",1000,0,2);
    hEffMap_JPsi  = new TH2F("hEffMap_JPsi","J/#psi counts (proxy for eff);p_{T};y",1000,0,10,1000,-3,3);
    hEffMap_D0    = new TH2F("hEffMap_D0","D^{0} counts (proxy for eff);p_{T};y",1000,0,10,1000,-3,3);
    // v2 profiles (unscaled)
    hV2JPsi = new TProfile("hV2JPsi","J/#psi #LTcos2#GT vs p_{T};p_{T};#LTcos2#GT",50,0,10);
    hV2D0   = new TProfile("hV2D0","D^{0} #LTcos2#GT vs p_{T};p_{T};#LTcos2#GT",50,0,10);
    // electron, kaion, poin like/unlike-sign spectra
    hJPsiMassVsPt1= new TH2F("hJPsiMassVsPt1","e^{+}e^{+} mass vs p_{T};M;p_{T}",1000,0,4,1000,0,10);
    hJPsiMassVsPt2= new TH2F("hJPsiMassVsPt2","e^{-}e^{-} mass vs p_{T};M;p_{T}",1000,0,4,1000,0,10);
    hJPsiMassVsPtULS  = new TH2F("hJPsiMassVsPtULS","e^{+}e^{-} mass vs p_{T};M;p_{T}",1000,0,4,1000,0,10);
    hD0MassVsPt1= new TH2F("hD0MassVsPt1","K^{+}#pi^{+} mass vs p_{T};M;p_{T}",1000,0,4,1000,0,10);
    hD0MassVsPt2= new TH2F("hD0MassVsPt2","K^{-}#pi^{-} mass vs p_{T};M;p_{T}",1000,0,4,1000,0,10);
    hD0MassVsPtULS  = new TH2F("hD0MassVsPtULS","K^{+}#pi^{-} mass vs p_{T};M;p_{T}",1000,0,4,1000,0,10);   
    //Event QA 
    hRefMultVz    = new TH2F("hRefMultVz","gRefMult vs Vz;Vz (cm);gRefMult",1000,-60,60,100,0,1000);
    mPsi2 = 0.f;
    mEpFinder = new StEpdEpFinder(10);
    mEpFinder->SetEpdHitFormat(2);
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

float StHFAnalysisMaker::trackBeta(const StPicoTrack* trk) const{
    constexpr float c_light=2.99792458e10; // cm/s
    int idx=trk->bTofPidTraitsIndex();
    float beta=std::numeric_limits<float>::quiet_NaN();
    auto picoDst = mPicoDstMaker->picoDst();
    if(idx>=0){
        const StPicoBTofPidTraits* tofPid = picoDst->btofPidTraits(idx);
        if(tofPid){
            beta = tofPid->btofBeta();
            if(beta<1e-4){
                TVector3 vtx3 = picoDst->event()->primaryVertex();
                StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
                TVector3 hit3 = tofPid->btofHitPos();
                StThreeVectorF hit(hit3.x(),hit3.y(),hit3.z());
                StPicoPhysicalHelix helix = trk->helix(picoDst->event()->bField());
                float L = tofPathLength(&vtx,&hit,helix.curvature());
                float tof = tofPid->btof();
                if(tof>0) beta = L/(tof*(c_light*1e-9));
            }
        }
    }
    return beta;
}

// --- D0 pairing ---
void StHFAnalysisMaker::runD0(){
    struct KPTrack{ const StPicoTrack* k; const StPicoTrack* p; TLorentzVector lv;};
    const double mK=0.493677, mPi=0.13957;
    auto makeLV=[&](const StPicoTrack* tr,double m){TVector3 p=tr->pMom(); return TLorentzVector(p, std::sqrt(p.Mag2()+m*m));};
    std::vector<const StPicoTrack*> kPlus=mKplus, kMinus=mKminus, piPlus=mPiplus, piMinus=mPiminus;

    auto pairLoop=[&](const std::vector<const StPicoTrack*>& K, const std::vector<const StPicoTrack*>& P,
                      TH1* hM, TH2* hMPt, bool sameCharge){
        for(size_t i=0;i<K.size();++i){
            const StPicoTrack* ka = K[i]; TLorentzVector lvK = makeLV(ka,mK);
            size_t jStart = sameCharge ? i+1 : 0;
            for(size_t j=jStart;j<P.size();++j){
                const StPicoTrack* pi = P[j];
                if(ka==pi) continue; // avoid identical track when same list
                TLorentzVector lvPi = makeLV(pi,mPi);
                TLorentzVector pair = lvK + lvPi;
                double mass=pair.M(); double pt=pair.Pt();
                if(hM)   hM->Fill(mass);
                if(hMPt) hMPt->Fill(mass,pt);
                if(!sameCharge){
                    double y = pair.Rapidity();
                    if(hD0PtY) hD0PtY->Fill(pt,y);
                    if(hEffMap_D0) hEffMap_D0->Fill(pt,y);
                    if(hPhiVsEP_D0 && hV2D0){
                        double phi = pair.Phi();
                        double dphi = TVector2::Phi_mpi_pi(phi - mPsi2);
                        hPhiVsEP_D0->Fill(pt,dphi);
                        hV2D0->Fill(pt,std::cos(2*dphi));
                    }
                }
            }
        }
    };

    pairLoop(kPlus, piPlus,  hD0BkgMass1, hD0MassVsPt1, true);
    pairLoop(kMinus,piMinus, hD0BkgMass2, hD0MassVsPt2, true);
    pairLoop(kPlus, piMinus, hD0Mass,     hD0MassVsPtULS, false);
}


void StHFAnalysisMaker::runHFE(){
    auto pico = mPicoDstMaker->picoDst();
    for(const auto* t : mElectrons){
        double p = t->pMom().Mag();
        int idx = t->bemcPidTraitsIndex();
        if(idx>=0){
            const StPicoBEmcPidTraits* bemc = pico->bemcPidTraits(idx);
            if(bemc){
                float e = bemc->bemcE();
                if(e>0){
                    float eop = e/p;
                    hEoverPvsP->Fill(p,eop);
                    hEOPInclusive->Fill(eop);
                    if(eop>HFCuts::PID::eopMin && eop<HFCuts::PID::eopMax) hNPEPt->Fill(p);
                }
            }
        }
    }
}

// --- Dielectron pairing ---
void StHFAnalysisMaker::runDielectronPairs(){
    struct ETrack{const StPicoTrack* tr; TLorentzVector lv;};
    std::vector<ETrack> plus, minus; plus.reserve(mElectrons.size()); minus.reserve(mElectrons.size());
    const double me=0.000511; // GeV
    for(const auto* t: mElectrons){
        TVector3 p = t->pMom();
        TLorentzVector lv(p, std::sqrt(p.Mag2()+me*me));
        ETrack et; et.tr=t; et.lv=lv;
        if(t->charge()>0) plus.push_back(et); else minus.push_back(et);
    }
    auto pairLoop=[&](const std::vector<ETrack>& A,const std::vector<ETrack>& B,
                      TH1* hM,TH2* hMPt,bool sameList){
        for(size_t i=0;i<A.size();++i){
            const TLorentzVector &l1 = A[i].lv;
            size_t jStart = sameList ? i+1 : 0;
            for(size_t j=jStart;j<B.size();++j){
                const TLorentzVector &l2 = B[j].lv;
                TLorentzVector pr = l1 + l2;
                if(hM)   hM->Fill(pr.M());
                if(hMPt) hMPt->Fill(pr.M(), pr.Pt());
                if(!sameList){
                    double pt = pr.Pt();
                    double y  = pr.Rapidity();
                    if(hJPsiPtY) hJPsiPtY->Fill(pt,y);
                    if(hPhiVsEP_JPsi && hV2JPsi){
                        double dphi = TVector2::Phi_mpi_pi(pr.Phi()-mPsi2);
                        hPhiVsEP_JPsi->Fill(pt,dphi);
                        hV2JPsi->Fill(pt, std::cos(2*dphi));
                                        }
                }
            }
        }
    };
    pairLoop(plus , plus , hJPsiBkgMass1 , hJPsiMassVsPt1 , true);
    pairLoop(minus, minus, hJPsiBkgMass2 , hJPsiMassVsPt2 , true);
    pairLoop(plus , minus, hJPsiMass   , hJPsiMassVsPtULS  , false);
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
        // TOF beta cuts
        auto pidOk=[&](const StPicoTrack* tr,float mass,float diffBeta,float diffMsq){
            float beta=trackBeta(tr);
            if(!(beta==beta)) return true; // no TOF
            float p=tr->pMom().Mag(); float betaExpected = p/std::sqrt(p*p+mass*mass);
            if(std::fabs(1.0/beta - 1.0/betaExpected) >= diffBeta) return false;
            float m2 = p*p*(1.0/(beta*beta)-1.0);
            return std::fabs(m2 - mass*mass) < diffMsq;
        };
        if(nSigE < HFCuts::PID::nSigmaE && pidOk(t,0.000511,HFCuts::PID::betaDiffPiMax,HFCuts::PID::msqDiffE)) mElectrons.push_back(t);
        if(nSigK < HFCuts::PID::nSigmaK && pidOk(t,0.493677,HFCuts::PID::betaDiffKMax,HFCuts::PID::msqDiffK)) ((t->charge()>0)? mKplus : mKminus).push_back(t);
        if(nSigPi< HFCuts::PID::nSigmaPi && pidOk(t,0.13957 ,HFCuts::PID::betaDiffPiMax,HFCuts::PID::msqDiffPi)) ((t->charge()>0)? mPiplus: mPiminus).push_back(t);
    }
    // event-level QA
    auto evt = picoEvt->event();
    if(evt) hRefMultVz->Fill(evt->primaryVertex().z(),evt->grefMult());

    // compute EP for this event
    const int kEpdHit=8;
    TClonesArray* epdHits = StPicoDst::picoArray(kEpdHit);
    StEpdEpInfo epInfo = mEpFinder->Results(epdHits, picoEvt->event()->primaryVertex(),0);
    mPsi2 = epInfo.FullPhiWeightedAndShiftedPsi(2);

    runD0(); runHFE(); runDielectronPairs();
    return kStOK;
}

Int_t StHFAnalysisMaker::Finish(){
    if(mOutFile.empty()){LOG_ERROR<<"Output file name is empty"<<endm; return kStFatal;}
    TFile *f = new TFile(mOutFile.c_str(),"RECREATE","HFAnalysisOutput",9);
    if(!f||f->IsZombie()){LOG_ERROR<<"Cannot create output file "<<mOutFile<<endm; return kStFatal;}

    

    // Let EPD finder write its own correction file (it changes gFile)
    if(mEpFinder) mEpFinder->Finish();

    // Ensure we are back in our output file directory
    f->cd();

    TH1* h1s[] = {
        hJPsiMass,
        hJPsiBkgMass1,
        hJPsiBkgMass2,
        hD0Mass,
        hD0BkgMass1,
        hD0BkgMass2,
        hNPEPt,
        hEOPInclusive,
        hED0_DeltaPhi
    };
    TH2* h2s[] = {
        hJPsiMassVsPt1,
        hJPsiMassVsPt2,
        hJPsiMassVsPtULS,
        hJPsiPtY,
        hPhiVsEP_JPsi,
        hEffMap_JPsi,
        hD0PtY,
        hEffMap_D0,
        hD0MassVsPt1,
        hD0MassVsPt2,
        hD0MassVsPtULS,
        hEoverPvsP,
        hPhiVsEP_D0,
        hRefMultVz
    };
    for(auto h:h1s) if(h) h->Write("",TObject::kOverwrite);
    for(auto h:h2s) if(h) h->Write("",TObject::kOverwrite);
    if(hV2JPsi) hV2JPsi->Write();
    if(hV2D0)   hV2D0->Write();

    f->Write();
    f->Close();
    return kStOK;
}
