// Macro to run the Pico HF analysis
#include <TSystem.h>
#include "StChain/StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StHFAnalysisMaker/StHFAnalysisMaker.h"

void runPicoHF(const char* inList="pico.list",
               const char* outFile="hfOut.root",
               int nEvents = 1000000000)
{
    gSystem->Load("libPicoHFAnalysis.so");

    StChain* chain = new StChain();
    StPicoDstMaker* pico = new StPicoDstMaker(0, inList, "picoDst");

    StHFAnalysisMaker* hf = new StHFAnalysisMaker("HF");
    hf->SetOutFile(outFile);
    hf->SetPicoDstMaker(pico);

    // Enable physics channels
    

    chain->Init();
    int iEvent = 0;
    while (chain->Make() == kStOk && iEvent < nEvents) {
        ++iEvent;
        if (iEvent % 10000 == 0) {
            std::cout << "Processed " << iEvent << " events" << std::endl;
        }
    }
    chain->Finish();
}
