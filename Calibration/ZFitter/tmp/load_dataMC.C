{
  gROOT->ProcessLine(".L macro/PlotDataMC.C+");
  gROOT->ProcessLine(".L macro/addChainWithFriends.C+");
  gROOT->ProcessLine(".L src/setTDRStyle.cc");
  setTDRStyle();

  TChain *data   = (TChain *) _file0->Get("selected");
  TChain *signalA = (TChain *) _file1->Get("selected");

  ReassociateFriends(_file0, data);
  ReassociateFriends(_file1, signalA);

   TDirectory *curDir = new TDirectory("curDir","");
   curDir->cd();

  TString outputPath="tmp/";

// look at standardDataMC.C for some examples  
}

