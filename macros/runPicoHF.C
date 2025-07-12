#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;

StChain *chain;

void runPicoHF(const char *list="pico.list",
               const char *out="hfOut.root",
               Long64_t    nEv=1e9)
{
gROOT->Macro("loadMuDst.C");
gSystem->Load("StDbBroker");
gSystem->Load("St_db_Maker.so");
gSystem->Load("StEmcUtil");
gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
loadSharedLibraries();

gSystem->Load("StPicoEvent");
gSystem->Load("StPicoDstMaker");

  gSystem->Load("StHFAnalysisMaker");                                 // your cons-built lib

  chain = new StChain("hfChain");
  StPicoDstMaker* pico = new StPicoDstMaker(2,list,"picoDstMaker");
//  StPicoDstMaker* pico = new StPicoDstMaker(2, list, "picoDstMaker");
  StHFAnalysisMaker* hf   = new StHFAnalysisMaker("HFMaker");
  hf->SetOutFile(out);
  hf->SetPicoDstMaker(pico);

  chain->Init();
  int nEntries = pico->chain()->GetEntries();
  cout<<"Processing "<<nEntries<<" events..."<<endl;
  for (int iEvent = 0; iEvent < nEntries; ++iEvent)
  {
    chain->Clear();
    if(iEvent && iEvent%10 == 0) cout<<"... finished processing "<<iEvent<<" events."<<endl;

    int iret = chain->Make();
    if (iret)
    {
      cout << "Bad return code!" << iret << endl;
      break;
    }
  }
  cout<<"Finished processing "<<nEntries<<" events."<<endl;

  chain->Finish();
  delete chain;

}