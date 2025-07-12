#ifndef STHFANALYSISMAKER_H
#define STHFANALYSISMAKER_H

#include "StMaker.h"
#include <string>
#include <vector>
#include "StEpdUtil/StEpdEpFinder.h"
class StPicoDstMaker;
class StPicoTrack;
class TH1F; class TH2F; class TProfile;

class StHFAnalysisMaker : public StMaker {
public:
    StHFAnalysisMaker(const char* name="HFMaker");
    virtual ~StHFAnalysisMaker();

    void    SetOutFile(const char* fn) { mOutFile = fn; }
    void    SetPicoDstMaker(StPicoDstMaker* mk) { mPicoDstMaker = mk; }

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

private:

    bool    passEventCuts();
    bool    goodTrack(const StPicoTrack* t);
    float   trackBeta(const StPicoTrack* trk) const;

    void    runJPsi();
    void    runD0();
    void    runHFE();
    void    runDielectronPairs();
    void    fitMassPeaks();
    void    SetPsi2(float psi){ mPsi2 = psi; }

    StPicoDstMaker* mPicoDstMaker = nullptr;
    std::string     mOutFile = "hfOut.root";

    TH1F* hJPsiMass = nullptr; TH2F* hJPsiPtY = nullptr;
    TH1F* hJPsiBkgMass1 = nullptr; TH1F* hJPsiBkgMass2 = nullptr; 
    TH2F* hJPsiMassVsPt1=nullptr,*hJPsiMassVsPt2=nullptr,*hJPsiMassVsPtULS=nullptr;
    TH1F* hD0Mass   = nullptr; TH2F* hD0PtY  = nullptr;
    TH1F* hD0BkgMass1 = nullptr; TH1F* hD0BkgMass2 = nullptr;
    TH1F* hNPEPt    = nullptr; TH2F* hEoverPvsP = nullptr;
    TH2F *hD0MassVsPt1=nullptr,*hD0MassVsPt2=nullptr,*hD0MassVsPtULS=nullptr;

    TH2F* hPhiVsEP_JPsi = nullptr;
    TH2F* hPhiVsEP_D0   = nullptr;
    TH1F* hED0_DeltaPhi = nullptr;
    TH1F* hEOPInclusive = nullptr;
    TH2F* hEffMap_JPsi = nullptr;
    TH2F* hEffMap_D0   = nullptr;
    // v2 profiles
    TProfile* hV2JPsi = nullptr;
    TProfile* hV2D0   = nullptr;
    
    std::vector<const StPicoTrack*> mElectrons, mKplus, mKminus, mPiplus, mPiminus;
    float mPsi2; 
    StEpdEpFinder* mEpFinder;

    TH2F* hRefMultVz   = nullptr;

    ClassDef(StHFAnalysisMaker,1)
};

#endif
