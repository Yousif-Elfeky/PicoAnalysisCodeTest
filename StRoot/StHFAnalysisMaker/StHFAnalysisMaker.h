#ifndef STHFANALYSISMAKER_H
#define STHFANALYSISMAKER_H

#include "StMaker.h"
#include <string>
#include <vector>
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
    // ===== internal helpers =====
    bool    passEventCuts();
    bool    goodTrack(const StPicoTrack* t);

    void    runJPsi();
    void    runD0();
    void    runHFE();
    void    runDielectronPairs();
    void    fitMassPeaks();
    void    SetPsi2(float psi){ mPsi2 = psi; }

    // ===== members =====
    StPicoDstMaker* mPicoDstMaker = nullptr;
    std::string     mOutFile = "hfOut.root";

    // JPsi
    TH1F* hJPsiMass = nullptr; TH2F* hJPsiPtY = nullptr;
    // D0
    TH1F* hD0Mass   = nullptr; TH2F* hD0PtY  = nullptr;
    // HFE
    TH1F* hNPEPt    = nullptr; TH2F* hEoverPvsP = nullptr;

    // ==== Additional observables ====
    // 3) Azimuthal anisotropy placeholders (phi-EP)
    TH2F* hPhiVsEP_JPsi = nullptr;
    TH2F* hPhiVsEP_D0   = nullptr;
    // 5) Correlation: e - D0 delta phi
    TH1F* hED0_DeltaPhi = nullptr;
    // 7) HFE E/p already covered; add inclusive electron cocktail QA
    TH1F* hEOPInclusive = nullptr;
    // 8) Acceptance / efficiency maps (raw counts here)
    TH2F* hEffMap_JPsi = nullptr;
    TH2F* hEffMap_D0   = nullptr;
    // Dielectron like/unlike-sign spectra
    TH1F *hMee_LSneg=nullptr,*hMee_LSpos=nullptr,*hMee_ULS=nullptr;
    // Fit functions
    TF1 *fJPsiSig=nullptr,*fJPsiBkg=nullptr,*fD0Sig=nullptr,*fD0Bkg=nullptr;
    TH2F *hMeePt_LSneg=nullptr,*hMeePt_LSpos=nullptr,*hMeePt_ULS=nullptr;
    // 9) Cached track lists per event to avoid rescanning
    std::vector<const StPicoTrack*> mElectrons, mKplus, mKminus, mPiplus, mPiminus;
    float mPsi2; ///< cached second-harmonic event-plane angle from EPD = 0.f;

    // Event QA-level: refMult vs Vz
    TH2F* hRefMultVz   = nullptr;

    ClassDef(StHFAnalysisMaker,1)
};

#endif
