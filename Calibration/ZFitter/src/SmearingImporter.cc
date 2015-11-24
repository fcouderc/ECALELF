#include "../interface/SmearingImporter.hh"
#include "../interface/BW_CB_pdf_class.hh"
#include <fstream>
#include <TTreeFormula.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TDirectory.h>
//#define DEBUG
#include <TStopwatch.h>
#define SELECTOR
#define FIXEDSMEARINGS
SmearingImporter::SmearingImporter(std::vector<TString> regionList, TString energyBranchName, TString commonCut):
  //  _chain(chain),
  _regionList(regionList),
  _scaleToy(1.01),
  _constTermToy(0.01),
  _energyBranchName(energyBranchName),
  _commonCut(commonCut),
  _eleID("loose"),
  _usePUweight(true),
  _useMCweight(true),
  _useR9weight(false),
  _usePtweight(false),
  _useZPtweight(false),
  _useFSRweight(false),
  _useWEAKweight(false),
  _excludeByWeight(true),
  _onlyDiagonal(false),
  _isSmearingEt(false),
  _pdfWeightIndex(0),
  cutter(false)
{
  cutter.energyBranchName=energyBranchName;
}


SmearingImporter::regions_cache_t SmearingImporter::GetCacheToy(Long64_t nEvents, bool isMC){
  std::cout<<"You are calling SmearingImporter::GetCacheToy. I didn't fix it with the general targetVariable. Be Careful"<<std::endl;
  regions_cache_t cache;
  TStopwatch myClock;
  myClock.Start();

  for(std::vector<TString>::const_iterator region_ele1_itr = _regionList.begin();
      region_ele1_itr != _regionList.end();
      region_ele1_itr++){
    for(std::vector<TString>::const_iterator region_ele2_itr = region_ele1_itr;
	region_ele2_itr != _regionList.end();
	region_ele2_itr++){
      event_cache_t eventCache;
      if(!(_onlyDiagonal && region_ele2_itr != region_ele1_itr))
	ImportToy(nEvents, eventCache, isMC);
      cache.push_back(eventCache);
    }
  }
  myClock.Stop();
  myClock.Print();
  return cache;
}

void SmearingImporter::ImportToy(Long64_t nEvents, event_cache_t& eventCache, bool isMC){
  std::cout<<"You are calling SmearingImporter::ImportToy. I didn't fix it with the general targetVariable. Be Careful"<<std::endl;
  TRandom3 gen(12345);
  if(isMC) nEvents*=2;
  //std::cout << "[DEBUG] constTermToy = " << _constTermToy << std::endl;
  for(Long64_t iEvent=0; iEvent < nEvents; iEvent++){
    ZeeEvent event;
    event.weight=1;
    event.energy_ele1=45;
    event.energy_ele2=45;
    //event.invMass = gen.BreitWigner(91.188,2.48);

    //event.targetVariable = gen.BreitWigner(91.188,2.48);
    //if(isMC==false) event.targetVariable*= sqrt(gen.Gaus(_scaleToy, _constTermToy) * gen.Gaus(_scaleToy, _constTermToy));//gen.Gaus(_scaleToy, _constTermToy); //
    eventCache.push_back(event);
  }
  return;
}



void SmearingImporter::Import(TTree *chain, regions_cache_t& cache, TString oddString, bool isMC, Long64_t nEvents, bool isToy, bool externToy){
  std::cout<<"[STATUS] Importing the events in SmearingImporter::Import"<<std::endl;

  std::ofstream ofs   ("global.txt", std::ofstream::app);
  std::ofstream ofs_0 ("BB.txt", std::ofstream::app);//BB
  std::ofstream ofs_1 ("BE.txt", std::ofstream::app);//BE
  std::ofstream ofs_2 ("EE.txt", std::ofstream::app);//EE

  TRandom3 gen(0);
  if(!isMC) gen.SetSeed(12345);
  TRandom3 excludeGen(12345);
  Long64_t excludedByWeight=0, includedEvents=0;

  // for the energy calculation
  Float_t         energyEle[2];
  Float_t         corrEle_[2]={1,1};
  Float_t         smearEle_[2]={1,1};
  bool hasSmearEle=false;

  TLorentzVector Ele1,Ele2;
  // for the angle calculation
  Float_t         etaEle[2];
  Float_t         phiEle[2];

  // for the weight calculation
  Float_t         weight=1.;
  Float_t         r9weight[2]={1,1};
  Float_t         ptweight[2]={1,1};
  Float_t         FSRweight=1.;
  Float_t         WEAKweight=1.;
  Float_t         zptweight[45]={1};
  Float_t         mcGenWeight=1;
  std::vector<double> *pdfWeights = NULL;

  Int_t           smearerCat[2];
  bool hasSmearerCat=false;

  // for toy repartition
  ULong64_t eventNumber;

  //
  Float_t ZPta;//check this

  //------------------------------
  chain->SetBranchAddress("eventNumber", &eventNumber);
  chain->SetBranchAddress("etaEle", etaEle);
  chain->SetBranchAddress("phiEle", phiEle);

  chain->SetBranchAddress(_energyBranchName, energyEle);
  std::cout<<"_energyBranchName in SmearingImporter is "<<_energyBranchName<<std::endl;
  if(chain->GetBranch("scaleEle")!=NULL){
    if(isToy==false || (externToy==true && isToy==true && isMC==false)){
    std::cout << "[STATUS] Adding electron energy correction branch from friend" << std::endl;
    chain->SetBranchAddress("scaleEle", corrEle_);
    cutter._corrEle=true;
    }
  } 

  if(chain->GetBranch("smearEle")!=NULL){
    if(isToy==false || (externToy==true && isToy==true && isMC==false)){
      std::cout << "[STATUS] Adding electron energy smearing branch from friend" << std::endl;
      if(isMC) chain->SetBranchAddress("smearSigmaEle", smearEle_);//e' la sigma della gaussiana
      else chain->SetBranchAddress("smearEle", smearEle_);//e' il numero gaussiano
      hasSmearEle=true;
    } 
  }

  if(chain->GetBranch("ZPta_energySCEle_regrCorrSemiParV5_ele")!=NULL){
    chain->SetBranchAddress("ZPta_energySCEle_regrCorrSemiParV5_ele",&ZPta);
  }
  if(!isMC && chain->GetBranch("pdfWeights_cteq66")!=NULL && _pdfWeightIndex>0){
    std::cout << "[STATUS] Adding pdfWeight_ctec66 branch from friend" << std::endl;
    chain->SetBranchAddress("pdfWeights_cteq66", &pdfWeights);
  }

  // the second term is to ensure that in case of toy study it is applied only to pseudo-data otherwise to MC
  // the third is only temporary because old ntuples does not have the branch
  // probably it will be needed in the future if the pdfSystematic branches are put in a separate tree
  if(_useFSRweight &&  isMC==false && (isToy==false || (externToy==true && isToy==true && isMC==false)) && chain->GetBranch("fsrWeight")!=NULL){
    std::cout << "[STATUS] Getting fsrWeight branch for tree: " << chain->GetTitle() << std::endl;
    chain->SetBranchAddress("fsrWeight", &FSRweight);
  }
  if(_useWEAKweight  && isMC==false && (isToy==false || (externToy==true && isToy==true && isMC==false)) && chain->GetBranch("weakWeight")!=NULL){
    std::cout << "[STATUS] Getting weakWeight branch for tree: " << chain->GetTitle() << std::endl;
    chain->SetBranchAddress("weakWeight", &WEAKweight);
  }

  if(chain->GetBranch("puWeight")!=NULL){
    std::cout << "[STATUS] Getting puWeight branch for tree: " << chain->GetTitle() << std::endl;
    chain->SetBranchAddress("puWeight", &weight);
  }

  if(chain->GetBranch("r9Weight")!=NULL){
    std::cout << "[STATUS] Getting r9Weight branch for tree: " << chain->GetTitle() << std::endl;
    chain->SetBranchAddress("r9Weight", r9weight);
  }  

  if(chain->GetBranch("ptWeight")!=NULL){
    std::cout << "[STATUS] Getting ptWeight branch for tree: " <<  chain->GetTitle() << std::endl;
    chain->SetBranchAddress("ptWeight", ptweight);
  }

  if(_useZPtweight && chain->GetBranch("ZPtWeight")!=NULL){
    std::cout << "[STATUS] Getting ZptWeight branch for tree: " <<  chain->GetTitle() << std::endl;
    chain->SetBranchAddress("ZPtWeight", zptweight);
  }

  if(chain->GetBranch("mcGenWeight")!=NULL){
    std::cout << "[STATUS] Getting mcGenWeight branch for tree: " <<  chain->GetTitle() << std::endl;
    chain->SetBranchAddress("mcGenWeight", &mcGenWeight);
  }

  if(chain->GetBranch("smearerCat")!=NULL){
    std::cout << "[STATUS] Getting smearerCat branch for tree: " <<  chain->GetTitle() << std::endl;
    chain->SetBranchAddress("smearerCat", smearerCat);
    hasSmearerCat=true;
  }

  if(hasSmearerCat==false){
    std::cerr << "[ERROR] Must have smearerCat branch" << std::endl;
    exit(1);
  }
  Long64_t entries = chain->GetEntryList()->GetN();
  std::cout<<"[DEBUG] GetEntryList()->GetN() del chain in Import "<<chain->GetEntryList()->GetN()<<std::endl;
  if(nEvents>0 && nEvents<entries){
    std::cout << "[INFO] Importing only " << nEvents << " events" << std::endl;
    entries=nEvents;
  }
  chain->LoadTree(chain->GetEntryNumber(0));
  Long64_t treenumber=-1;


  std::vector< std::pair<TTreeFormula*, TTreeFormula*> > catSelectors;
  if(hasSmearerCat==false){
    for(std::vector<TString>::const_iterator region_ele1_itr = _regionList.begin();
	region_ele1_itr != _regionList.end();
	region_ele1_itr++){
      for(std::vector<TString>::const_iterator region_ele2_itr = region_ele1_itr;
	  region_ele2_itr != _regionList.end();
	  region_ele2_itr++){
	
	if(region_ele2_itr==region_ele1_itr){
	  TString region=*region_ele1_itr;
	  region.ReplaceAll(_commonCut,"");
	  TTreeFormula *selector = new TTreeFormula("selector"+(region), cutter.GetCut(region+oddString, isMC), chain);
	  catSelectors.push_back(std::pair<TTreeFormula*, TTreeFormula*>(selector,NULL));
	  //selector->Print();
	} else if(!_onlyDiagonal){
	  TString region1=*region_ele1_itr;
	  TString region2=*region_ele2_itr;
	  region1.ReplaceAll(_commonCut,"");
	  region2.ReplaceAll(_commonCut,"");
	  TTreeFormula *selector1 = new TTreeFormula("selector1"+region1+region2, 
						     cutter.GetCut(region1+oddString, isMC, 1) && 
						     cutter.GetCut(region2+oddString, isMC, 2),
						     chain);
	  TTreeFormula *selector2 = new TTreeFormula("selector2"+region1+region2,
						     cutter.GetCut(region1+oddString, isMC, 2) && 
						     cutter.GetCut(region2+oddString, isMC, 1),
						     chain);
	  catSelectors.push_back(std::pair<TTreeFormula*, TTreeFormula*>(selector1,selector2));
	  //selector1->Print();
	  //selector2->Print();
	  
	} else catSelectors.push_back(std::pair<TTreeFormula*, TTreeFormula*>(NULL, NULL));
	
      }
    }
  }//hasSmearerCat false

  TRandom3 rand(0);//Use a seed generated using machine clock (different every second)
  double rn;
  int test_entries=5000;
  std::cout<<"Starting the loop in Import"<<std::endl;
  std::cout<<"entries are "<<entries<<std::endl;
  for(Long64_t jentry=0; jentry < entries; jentry++){
    //for(Long64_t jentry=0; jentry < test_entries; jentry++){//std::cout<<"WARNING!! you are running on test_entries in SmearingImporter::Import()"<<std::endl;
    //std::cout<<"jentry is "<<jentry<<std::endl;
    Long64_t entryNumber= chain->GetEntryNumber(jentry);
    chain->GetEntry(entryNumber);
    if(isToy){
      std::cout<<"insideToys"<<std::endl;
      int divider=20;
      int modulo=eventNumber%divider;
      if(jentry<10){
	std::cout<<"[toyMC INFO] In SmearingImporter::Import samples splitted % "<<divider<<std::endl;
	std::cout << "Dividing toyMC events: " << isMC << "\t" << eventNumber << "\t" << modulo
		  << "\t" << mcGenWeight 
		  << std::endl;
	
      }
      int data_index=18;
      if(data_index > 1000){
	//closure tests
	//data --> 0,1 /5 //perche' c'e' continue
	//MC   --> 2,3,4/5
	if(isMC && modulo<2) continue;
	if(!isMC && modulo>=2) continue;
      }else{
	//Stat_check study
	//data are identified by this data_index, the rest is MC
	if(isMC && modulo==data_index) continue; //MC
	if(!isMC && modulo!=data_index) continue; //dati
      }
    }//isToy

    //std::cout<<"in SmearingImporter after toys"<<std::endl;
    // reject events:
    if(weight>3) continue;

    if (hasSmearerCat==false && chain->GetTreeNumber() != treenumber) {
      treenumber = chain->GetTreeNumber();
      for(std::vector< std::pair<TTreeFormula*, TTreeFormula*> >::const_iterator catSelector_itr = catSelectors.begin();
	  catSelector_itr != catSelectors.end();
	  catSelector_itr++){
	
	catSelector_itr->first->UpdateFormulaLeaves();
	if(catSelector_itr->second!=NULL)       catSelector_itr->second->UpdateFormulaLeaves();
      }
    }

    int evIndex=-1;
    bool _swap=false;
    //std::cout<<"evIndex is "<<evIndex<<std::endl;
    if(!hasSmearerCat){
      for(std::vector< std::pair<TTreeFormula*, TTreeFormula*> >::const_iterator catSelector_itr = catSelectors.begin();
	  catSelector_itr != catSelectors.end();
	  catSelector_itr++){
	_swap=false;
	TTreeFormula *sel1 = catSelector_itr->first;
	TTreeFormula *sel2 = catSelector_itr->second;
	if(sel1==NULL) continue; // is it possible?
	if(sel1->EvalInstance()==false){
	  if(sel2==NULL || sel2->EvalInstance()==false) continue;
	  else _swap=true;
	}
	evIndex=catSelector_itr-catSelectors.begin();
      }
    }else{
      evIndex=smearerCat[0];
      _swap=smearerCat[1];
      if(jentry<2) std::cout << evIndex << "\t" << _swap << std::endl;
    }
    //std::cout<<"evIndex is "<<evIndex<<std::endl;
    if(evIndex<0) continue; // event in no category
    
    ZeeEvent event;
    
    float t1=TMath::Exp(-etaEle[0]);
    float t2=TMath::Exp(-etaEle[1]);
    float t1q = t1*t1;
    float t2q = t2*t2;
    
    if(isMC && hasSmearEle){//MC is smeared with a Gaussian factor
      //std::cout<<"smearEle[0] is "<<smearEle_[0]<<std::endl;
      smearEle_[0]=gen.Gaus(1,smearEle_[0]);
      smearEle_[1]=gen.Gaus(1,smearEle_[1]);
      //std::cout<<"smearEle[0] after is "<<smearEle_[0]<<std::endl;
    }

    //------------------------------
    Ele1.SetPtEtaPhiE(event.energy_ele1/cosh(etaEle[0]),etaEle[0],phiEle[0],event.energy_ele1);
    Ele2.SetPtEtaPhiE(event.energy_ele2/cosh(etaEle[1]),etaEle[1],phiEle[1],event.energy_ele2);
    float ZPt=(Ele1+Ele2).Pt();
    
    if(_swap){
      event.energy_ele2 = energyEle[0] * corrEle_[0] * smearEle_[0];
      event.energy_ele1 = energyEle[1] * corrEle_[1] * smearEle_[1];
    } else {
      event.energy_ele1 = energyEle[0] * corrEle_[0] * smearEle_[0];
      event.energy_ele2 = energyEle[1] * corrEle_[1] * smearEle_[1];
    }
    
    Ele1.SetPtEtaPhiE(event.energy_ele1/cosh(etaEle[0]),etaEle[0],phiEle[0],event.energy_ele1);
    Ele2.SetPtEtaPhiE(event.energy_ele2/cosh(etaEle[1]),etaEle[1],phiEle[1],event.energy_ele2);
    float ZPt_after=(Ele1+Ele2).Pt();//corrections applied
    //std::cout<<"I am doing something ZPt is <<"<<ZPt_after<<std::endl;
    if(ZPt_after > 10) continue; //only events with Pt < 10 GeV

    std::ofstream smear_file_BB ("smear_applied_BB.txt", std::ofstream::app);
    //std::ofstream smear_file_BE ("smear_applied_BE.txt", std::ofstream::app);
    std::ofstream smear_file_EE ("smear_applied_EE.txt", std::ofstream::app);
    std::ofstream corr_file_BB ("corr_applied_BB.txt", std::ofstream::app);
    //std::ofstream corr_file_BE ("smear_applied_BE.txt", std::ofstream::app);
    std::ofstream corr_file_EE ("corr_applied_EE.txt", std::ofstream::app);
    if(evIndex==0){
      smear_file_BB << smearEle_[0]<<std::endl;
      smear_file_BB << smearEle_[1]<<std::endl;
      corr_file_BB << corrEle_[0]<<std::endl;
      corr_file_BB << corrEle_[1]<<std::endl;
    }else if(evIndex==2){
      smear_file_EE << smearEle_[0]<<std::endl;
      smear_file_EE << smearEle_[1]<<std::endl;
      corr_file_EE << corrEle_[0]<<std::endl;
      corr_file_EE << corrEle_[1]<<std::endl;
    }

    //Different pt scenarios
    //if((ZPt_after < 10) || (ZPt_after>20)) continue;
    //if((ZPt_after < 20) || (ZPt_after>40)) continue;
    //if((ZPt_after < 40) || (ZPt_after>100)) continue;
    //if((ZPt_after < 100) || (ZPt_after>200)) continue;

    //ZPta doesn't match ZPt_after!! Check this
    /*if(ZPt_after>10 && evIndex>0){//Useful for the future
      //check also event.energy_ele1
    std::cout<<"ZPt before correction is "<<ZPt<<std::endl;//invariant under swap
    std::cout<<"ZPt after correction is "<<ZPt_after<<std::endl;//invariant under swap
    std::cout<<"ZPta is "<<ZPta<<std::endl;//invariant under swap
    std::cout<<"evIndex is "<<evIndex<<std::endl;
    std::cout<<"corrEle[0] is "<<corrEle_[0]<<std::endl;
    std::cout<<"corrEle[1] is "<<corrEle_[1]<<std::endl;
    std::cout<<"smearEle[0] is "<<smearEle_[0]<<std::endl;
    std::cout<<"smearEle[1] is "<<smearEle_[0]<<std::endl;
    }*/
    //loop to fill targetVariable (it's a vector now)
    std::vector<TString> targets=GetTargetVariable();
    for(int var=0;((unsigned)var)<targets.size();var++){
      (event.targetVariable).push_back(0.);
    }

    for(int var=0;((unsigned)var)<targets.size();var++){
      event.targetVariable[var]=-999;
      if(targets[var]=="invMass"){
	// to calculate the invMass: invMass = sqrt(2 * energy_ele1 * energy_ele2 * angle_eta_ele1_ele2)
	//if(event.invMass < 70 || event.invMass > 110) continue;
	//event.angle_eta_ele1_ele2=  (1-((1-t1q)*(1-t2q)+4*t1*t2*cos(phiEle[0]-phiEle[1]))/((1+t1q)*(1+t2q)));
	event.targetVariable[var]= sqrt(2 * event.energy_ele1 * event.energy_ele2 *
					(1-((1-t1q)*(1-t2q)+4*t1*t2*cos(phiEle[0]-phiEle[1]))/((1+t1q)*(1+t2q)))
					);
      }else if(targets[var]=="pt2Sum"){
	double pt1=event.energy_ele1/cosh(etaEle[0]);
	double pt2=event.energy_ele2/cosh(etaEle[1]);
	
	if(GetConfiguration()=="leading"){
	  if(pt1 > pt2){
	    event.pt1=pt1;
	    event.pt2=pt2;
	    event.isSwapped=false;
	  }else{
	    event.pt1=pt2;
	    event.pt2=pt1;
	    event.isSwapped=true;
	  }
	}else if(GetConfiguration()=="random"){
	  rn=rand.Uniform(0,1);
	  if(rn<0.5){
	    event.pt1=pt1;
	    event.pt2=pt2;
	    event.isSwapped=false;
	  }else{
	    event.pt1=pt2;
	    event.pt2=pt1;
	    event.isSwapped=true;
	  }
	}
	
	event.targetVariable[var]=sqrt(pow(event.pt1,2)+pow(event.pt2,2));
	
	if(event.isSwapped){
	  ofs<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[1]<<" "<<etaEle[0]<<std::endl;
	}else{
	  ofs<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[0]<<" "<<etaEle[1]<<std::endl;
	}
	
	if(evIndex==0){
	  if(event.isSwapped){
	    ofs_0<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[1]<<" "<<etaEle[0]<<std::endl;
	  }else{
	    ofs_0<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[0]<<" "<<etaEle[1]<<std::endl;
	  }
	}else if(evIndex==1){
	  if(event.isSwapped){
	    ofs_1<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[1]<<" "<<etaEle[0]<<std::endl;
	  }else{
	    ofs_1<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[0]<<" "<<etaEle[1]<<std::endl;
	  }
	}else if(evIndex==2){
	  if(event.isSwapped){
	    ofs_2<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[1]<<" "<<etaEle[0]<<std::endl;
	  }else{
	    ofs_2<<event.targetVariable[var]<<" "<<event.pt1<<" "<<event.pt2<<" "<<etaEle[0]<<" "<<etaEle[1]<<std::endl;
	  }
	}

      }else if(targets[var]=="ptRatio"){
	//if(ZPt_after>10){
	//event.targetVariable[var]=-999;Checked: it works because Integral() avoids underflows and overflows
	//}else{
	
	//Decide which is 1 and which is 2, in this case
	//Option 1: 1 is the leading and 2 is the subleading
	double pt1=event.energy_ele1/cosh(etaEle[0]);
	double pt2=event.energy_ele2/cosh(etaEle[1]);
	
	if(GetConfiguration()=="leading"){
	  if(pt1 > pt2){
	    event.pt1=pt1;
	    event.pt2=pt2;
	    event.isSwapped=false;
	  }else{
	    event.pt1=pt2;
	    event.pt2=pt1;
	    event.isSwapped=true;
	  }
	}else if(GetConfiguration()=="random"){
	  rn=rand.Uniform(0,1);
	  if(rn<0.5){
	    event.pt1=pt1;
	    event.pt2=pt2;
	    event.isSwapped=false;
	  }else{
	    event.pt1=pt2;
	    event.pt2=pt1;
	    event.isSwapped=true;
	  }
	}
	
	event.targetVariable[var]=(event.pt1/event.pt2);
	ofs<<event.targetVariable[var]<<std::endl;
	if(evIndex==0){
	  ofs_0<<event.targetVariable[var]<<" ";
	}else if(evIndex==1){
	  ofs_1<<event.targetVariable[var]<<" ";
	}else if(evIndex==2){
	  ofs_2<<event.targetVariable[var]<<" ";
	}
      }else if(targets[var]=="ptSum"){
	double pt1=event.energy_ele1/cosh(etaEle[0]);
	double pt2=event.energy_ele2/cosh(etaEle[1]);
	if(GetConfiguration()=="leading"){
	  if(pt1 > pt2){
	    event.pt1=pt1;
	    event.pt2=pt2;
	    event.isSwapped=false;
	  }else{
	    event.pt1=pt2;
	    event.pt2=pt1;
	    event.isSwapped=true;
	  }
	}else if(GetConfiguration()=="random"){
	  rn=rand.Uniform(0,1);
	  if(rn<0.5){
	    event.pt1=pt1;
	    event.pt2=pt2;
	    event.isSwapped=false;
	  }else{
	    event.pt1=pt2;
	    event.pt2=pt1;
	    event.isSwapped=true;
	  }
	}

	event.targetVariable[var]=(event.pt1 + event.pt2);

	//event.targetVariable[var]=(pt1 + pt2);
	//event.pt1=pt1;//Needed for RooSmearer::SetSmearedHisto
	//event.pt2=pt2;
      }
    }//loop over targetVariables

    if(_isSmearingEt){
      if(_swap){
	event.energy_ele2/=cosh(etaEle[0]);
	event.energy_ele1/=cosh(etaEle[1]);
      } else{
	event.energy_ele1/=cosh(etaEle[0]);
	event.energy_ele2/=cosh(etaEle[1]);
      }	
    }

    event.weight = 1.;
    if(_usePUweight) event.weight *= weight;
    if(isMC){//remove this check
      if(weight==0){
	event.weight=1.; //puWeight seems to be 0
      }
    }    
    if(_useR9weight) event.weight *= r9weight[0]*r9weight[1];
    if(_usePtweight) event.weight *= ptweight[0]*ptweight[1];
    if(_useFSRweight) event.weight *= FSRweight;
    if(_useWEAKweight) event.weight *= WEAKweight;
    if(_useZPtweight && isMC && _pdfWeightIndex>0) event.weight *= zptweight[_pdfWeightIndex];

    if(!isMC && _pdfWeightIndex>0 && pdfWeights!=NULL){
      if(((unsigned int)_pdfWeightIndex) > pdfWeights->size()) continue;
      event.weight *= ((*pdfWeights)[0]<=0 || (*pdfWeights)[0]!=(*pdfWeights)[0] || (*pdfWeights)[_pdfWeightIndex]!=(*pdfWeights)[_pdfWeightIndex])? 0 : (*pdfWeights)[_pdfWeightIndex]/(*pdfWeights)[0];
    }else{
      if(!isMC && _pdfWeightIndex>0){
	std::cerr << "[ERROR] requested pdfWeights but not set by getentry" << std::endl;
	std::cerr << "[ERROR] jentry=" << jentry << "; chain: " << chain->GetName() << "\t" << chain->GetTitle()  << std::endl;
	if(jentry<10) continue;
	else exit(1);
      }
    }

    //As simple as that                                                                                                                                                        
    if(isMC){
      if(mcGenWeight > 0){
        event.weight *=1;
      }else{
        event.weight *=-1;
      }
    }

    //This was for pythia
//    if(mcGenWeight != -1){
//      if(_useMCweight && !_excludeByWeight) event.weight *= mcGenWeight;
//      
//      if(_excludeByWeight && mcGenWeight!=1){
//	float rnd = excludeGen.Rndm();
//	if(jentry < 10)  std::cout << "mcGen = " << mcGenWeight << "\t" << rnd << std::endl;
//	if(mcGenWeight < rnd){ /// \todo fix the reject by weight in case of pu,r9,pt reweighting
//	  
//	  excludedByWeight++;
//	  continue;
//	}// else event.weight=1;
//      }
//    }

#ifdef DEBUG      
    if(jentry<10 || event.weight!=event.weight || event.weight>2){
      std::cout << "jentry = " << jentry 
		<< "\tevent.weight = " << event.weight 
		<< "\t" << weight << "\t" << mcGenWeight
		<< "\t" << r9weight[0] << " " << r9weight[1] 
		<< "\t" << ptweight[0] << " " << ptweight[1]
		<< "\t" << zptweight[0] 
		<< "\t" << WEAKweight << "\t" << FSRweight
		<< std::endl;
    }
#endif
    //if(event.weight<=0 || event.weight!=event.weight || event.weight>10) continue;
    if(event.weight>10) continue;//negative weights are now possible

#ifdef FIXEDSMEARINGS
    if(isMC){
      event.smearings_ele1 = new float[NSMEARTOYLIM];
      event.smearings_ele2 = new float[NSMEARTOYLIM];
      for(int i=0; i < NSMEARTOYLIM; i++){
	event.smearings_ele1[i] = (float) gen.Gaus(0,1);
	event.smearings_ele2[i] = (float) gen.Gaus(0,1);
      }
    }else{
      event.smearings_ele1 = new float[1];
      event.smearings_ele2 = new float[1];
      event.smearings_ele1[0] = (float) gen.Gaus(0,1);
      event.smearings_ele2[0] = (float) gen.Gaus(0,1);
    }	
#endif
    includedEvents++;
    cache.at(evIndex).push_back(event);
    //(cache[evIndex]).push_back(event);
  }//loop over entries

  std::cout << "[INFO] Importing events: " << includedEvents << "; events excluded by weight: " << excludedByWeight << std::endl;
  chain->ResetBranchAddresses();
  chain->GetEntry(0);
  return;
}

SmearingImporter::regions_cache_t SmearingImporter::GetCache(TChain *_chain, bool isMC, bool odd, Long64_t nEvents, bool isToy, bool externToy){
  std::cout<<"[STATUS] Preparing the cache and activating branches in SmearingImporter::GetCache"<<std::endl;
  //std::cout << "[STATUS] --- Setting toy cache for data" << std::endl;
  //data_events_cache = importer.GetCache(_signal_chain, false, false, nEvents, true, externToy); //importer.GetCacheToy(nEvents, false);
  TString eleID_="eleID_"+_eleID;

  TString oddString;
  if(odd) oddString+="-odd";
  regions_cache_t cache;
  TStopwatch myClock;
  myClock.Start();
  _chain->GetEntry(0);

  _chain->SetBranchStatus("*", 0);

  _chain->SetBranchStatus("etaEle", 1);
  _chain->SetBranchStatus("phiEle", 1);
  _chain->SetBranchStatus(_energyBranchName, 1);
  if(isToy) _chain->SetBranchStatus("eventNumber",1);
  //  std::cout << _chain->GetBranchStatus("seedXSCEle") <<  std::endl;
  //  std::cout << _chain->GetBranchStatus("etaEle") <<  std::endl;

  if(_chain->GetBranch("scaleEle")!=NULL){
    std::cout << "[STATUS] Activating branch scaleEle" << std::endl;
    _chain->SetBranchStatus("scaleEle", 1);
    cutter._corrEle=true;
  }

  if(_chain->GetBranch("smearEle")!=NULL){
    std::cout << "[STATUS] Activating branch smearEle" << std::endl;
    _chain->SetBranchStatus("smearEle", 1);
  } 

 
  if(_chain->GetBranch("r9Weight")!=NULL)  _chain->SetBranchStatus("r9Weight", 1);
  if(_chain->GetBranch("puWeight")!=NULL)  _chain->SetBranchStatus("puWeight", 1);
  if(_chain->GetBranch("ptWeight")!=NULL)  _chain->SetBranchStatus("ptWeight", 1);
  if(_chain->GetBranch("mcGenWeight")!=NULL)  _chain->SetBranchStatus("mcGenWeight", 1);
  if(_chain->GetBranch("smearerCat")!=NULL){
    //std::cout << "[STATUS] Activating branch smearerCat" << std::endl;
    _chain->SetBranchStatus("smearerCat", 1);
  }


  for(std::vector<TString>::const_iterator region_ele1_itr = _regionList.begin();
      region_ele1_itr != _regionList.end();
      region_ele1_itr++){
    std::set<TString> branchNames = cutter.GetBranchNameNtuple(_commonCut+"-"+eleID_+"-"+*region_ele1_itr);
    for(std::set<TString>::const_iterator itr = branchNames.begin();
	itr != branchNames.end(); itr++){
      _chain->SetBranchStatus(*itr, "1");
    }
  }
 
  TString evListName="evList_";
  evListName+=_chain->GetTitle();
  evListName+="_all";
  TEntryList *oldList = _chain->GetEntryList();
  if(oldList==NULL){
    std::cout << "[STATUS] Setting entry list (commotCut screening): " << evListName << std::endl;
    std::cout << "[STATUS] Marking events if "<<cutter.GetCut(_commonCut+"-"+eleID_,isMC)<<std::endl;
    _chain->Draw(">>"+evListName, cutter.GetCut(_commonCut+"-"+eleID_,isMC), "entrylist");
    //_chain->Draw(">>"+evListName, "", "entrylist");
    TEntryList *elist_all = (TEntryList*)gDirectory->Get(evListName);
    //  elist_all->SetBit(!kCanDelete);
    _chain->SetEntryList(elist_all);
  }
  for(std::vector<TString>::const_iterator region_ele1_itr = _regionList.begin();
      region_ele1_itr != _regionList.end();
      region_ele1_itr++){
    for(std::vector<TString>::const_iterator region_ele2_itr = region_ele1_itr;
	region_ele2_itr != _regionList.end();
	region_ele2_itr++){
      
      event_cache_t eventCache;
      cache.push_back(eventCache);
    }
  }
  std::cout<<"nEvents just before Import (smearingImporeter::GetCache) is "<<nEvents<<std::endl;
  Import(_chain, cache, oddString, isMC, nEvents, isToy, externToy);
#ifdef DEBUG
  int index=0;
  for(std::vector<TString>::const_iterator region_ele1_itr = _regionList.begin();
      region_ele1_itr != _regionList.end();
      region_ele1_itr++){
    for(std::vector<TString>::const_iterator region_ele2_itr = region_ele1_itr;
	region_ele2_itr != _regionList.end();
	region_ele2_itr++){
      std::cout << "[INFO] Category " << index << " " <<  *region_ele1_itr << *region_ele2_itr
		<< " filled with " << cache[index].size() << " entries" 
		<< std::endl;
      index++;
    }
  }
#endif
  myClock.Stop();
  myClock.Print();
  return cache;
}


