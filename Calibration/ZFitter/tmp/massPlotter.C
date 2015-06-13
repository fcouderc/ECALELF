{
  TString OutputPath="tmp/Plots/";
  std::vector<TString> mcLabel_vec;
  mcLabel_vec.push_back("Madgraph");
  c = PlotDataMCs(data, MakeChainVector(signalA), "invMass_SC_regrCorrSemiParV5_ele", "(100,80,100)", "eleID_loose-trigger-noPF-EB", "", "DATA", mcLabel_vec, "Mass [GeV]","", false, true,false, false, false);
  c->SaveAs(OutputPath+"Mass_raw.png"); 

  //Smearing and Scale applied with (....,true,true)
  c = PlotDataMCs(data, MakeChainVector(signalA), "invMass_SC_regrCorrSemiParV5_ele", "(100,80,100)", "eleID_loose-trigger-noPF-EB", "", "DATA", mcLabel_vec, "Mass [GeV]","", false, true,false, true, true);
  c->SaveAs(OutputPath+"Mass_corrected.png");
}
