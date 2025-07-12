#ifndef STHFANALYSISMAKER_H
#define STHFANALYSISMAKER_H

#include "StMaker.h"
#include <string>
#include <vector>
#include "StEpdUtil/StEpdEpFinder.h"
class StPicoDstMaker;
class StPicoTrack;
class TH1F; class TH2F;

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

    void    runJPsi();
    void    runD0();
    void    runHFE();
    void    runDielectronPairs();
    void    fitMassPeaks();
    void    SetPsi2(float psi){ mPsi2 = psi; }

    StPicoDstMaker* mPicoDstMaker = nullptr;
    std::string     mOutFile = "hfOut.root";

    TH1F* hJPsiMass = nullptr; TH2F* hJPsiPtY = nullptr;
    TH1F* hD0Mass   = nullptr; TH2F* hD0PtY  = nullptr;
    TH1F* hNPEPt    = nullptr; TH2F* hEoverPvsP = nullptr;

    TH2F* hPhiVsEP_JPsi = nullptr;
    TH2F* hPhiVsEP_D0   = nullptr;
    TH1F* hED0_DeltaPhi = nullptr;
    TH1F* hEOPInclusive = nullptr;
    TH2F* hEffMap_JPsi = nullptr;
    TH2F* hEffMap_D0   = nullptr;
    TH1F *hMee_LSneg=nullptr,*hMee_LSpos=nullptr,*hMee_ULS=nullptr;
    TF1 *fJPsiSig=nullptr,*fJPsiBkg=nullptr,*fD0Sig=nullptr,*fD0Bkg=nullptr;
    TH2F *hMeePt_LSneg=nullptr,*hMeePt_LSpos=nullptr,*hMeePt_ULS=nullptr;
    std::vector<const StPicoTrack*> mElectrons, mKplus, mKminus, mPiplus, mPiminus;
    float mPsi2; 
    StEpdEpFinder* mEpFinder;

    TH2F* hRefMultVz   = nullptr;

    ClassDef(StHFAnalysisMaker,1)
};

#endif
