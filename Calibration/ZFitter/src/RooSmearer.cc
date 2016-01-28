#include "../interface/RooSmearer.hh"
//#include <RooPullVar.h>
#include <TSystem.h>
#include <TIterator.h>
#include <TFile.h>
#define FIXEDSMEARINGS
//#define DEBUG
//#define CPU_DEBUG
//#define FUNC_DEBUG
//#define MEM_DEBUG

#define debug_vector
RooSmearer::~RooSmearer(void){
  
}

RooSmearer::RooSmearer(const RooSmearer& old, const char* newname){
  //  RooSmearer(newname, old._data_chain, old._signal_chain, NULL,
  //	     old.importer._regionList, old._paramSet, old.importer._energyBranchName);
  RooSmearer();
}
  
RooSmearer::RooSmearer(const char *name,  ///< name of the variable
		       TChain *data_chain_,  
		       TChain *signal_chain_,
		       TChain *bkg_chain_, 
		       std::vector<TString> regionList, ///<list of the region names (used for the categorization)
		       std::vector<RooArgSet> params, ///< parameters divided by category (also rooformulavar)
		       RooArgSet parset, ///< list of real variables 
		       TString energyBranchName
		       ):
  RooAbsReal(name,energyBranchName),
  importer(regionList, energyBranchName),
  _params_vec(params),
  _paramSet("paramSet","Set of parameters",this),
   //targetVariable_min_(-999), targetVariable_max_(+999), targetVariable_bin_(2000),//default values for ranges, now those are vector<float>
  //targetVariable_("invMass"),//default target Variable is invMass, now is a vector<string>
  deltaNLLMaxSmearToy(330),
  _deactive_minEventsDiag(1000), _deactive_minEventsOffDiag(1500), _nSmearToy(20), 
  nllBase(0),
  nllVar("nll","",0,1e20),
  _isDataSmeared(false),
  _autoBin(false),
  _autoNsmear(false),
  smearscan(false),
  dataset(NULL)
{
  _data_chain = data_chain_;
  _signal_chain = signal_chain_;

  TIterator *it = parset.createIterator();
  for(RooRealVar *v = (RooRealVar*)it->Next(); v!=NULL; v = (RooRealVar*)it->Next()){
    if( v->isLValue()){ //! v->isConstant() &&) 
      _paramSet.add((*v));
    }
  }


  rgen_=new TRandom3(0); // inizializzo il generatore con seed impostato sull'ora
  myClock=new TStopwatch();
  
  _markov.SetParameters(_paramSet);
  nllMin=1e20;
}


void RooSmearer::SetCache(Long64_t nEvents, bool cacheToy, bool externToy){
  //regions_cache_t GetCache(TChain *_chain, bool isMC, bool odd, Long64_t nEvents=0, bool isToy=false, bool externToy=true, TString targetVariable="invMass");
  std::cout << "------------------------------------------------------------" << std::endl;
  std::cout << "[STATUS] Importing cache events in RooSmearer::SetCache" << std::endl;
  std::cout << "[STATUS] in RooSmearer::SetCache nEvents is " <<nEvents<< std::endl;

  if(data_events_cache.empty()){
    if(cacheToy){
      std::cout << "[STATUS] --- Setting toy cache for data in RooSmearer::SetCache" << std::endl;
      data_events_cache = importer.GetCache(_signal_chain, false, false, nEvents, true, externToy); //importer.GetCacheToy(nEvents, false);
    }else {
      std::cout << "[STATUS] --- Setting cache for data in RooSmearer::SetCache" << std::endl;
      data_events_cache = importer.GetCache(_data_chain, false, false, nEvents);
    }
  }
  if(!mc_events_cache.empty()){
    std::cerr << "[ERROR] mc_events_cache not empty: " << mc_events_cache.size() << std::endl;
    exit(1);
  }
  if(mc_events_cache.empty()){
    if(cacheToy){
      std::cout << "[STATUS] --- Setting toy cache for mc in RooSmearer::SetCache" << std::endl;
      mc_events_cache = importer.GetCache(_signal_chain, true, false, nEvents,true, externToy); //importer.GetCacheToy(nEvents, true);
    }else {
      std::cout << "[STATUS] --- Setting cache for mc in RooSmearer::SetCache" << std::endl;
      mc_events_cache = importer.GetCache(_signal_chain, true, false, nEvents);
    }
  }


  //#ifdef DEBUG
    std::cout << "[DEBUG] RooSmearer::SetCache --> Data events size: " << data_events_cache.size() << std::endl;
    std::cout << "[DEBUG] RooSmearer::SetCache --> MC events size: " << data_events_cache.size() << std::endl;
    //#endif

  return;
}

void RooSmearer::InitCategories(bool mcToy){ //bool mcToy not used before
  std::cout<<"[STATUS] In RooSmearer::InitCategories "<<std::endl;
  int index=0;
  //Memo: (int)3=2 : if you define 2 cats on single ele, you'll have (n)*(n+1)/2 "combined cats"; "+1 is meant to castrate (int)"
  ZeeCategories.reserve((int)(importer._regionList.size()*(importer._regionList.size()+1)/2 +1));
  for(std::vector<TString>::const_iterator region_ele1_itr = importer._regionList.begin();
      region_ele1_itr != importer._regionList.end();
      region_ele1_itr++){
    for(std::vector<TString>::const_iterator region_ele2_itr = region_ele1_itr; //importer._regionList.begin();
	region_ele2_itr != importer._regionList.end();
	region_ele2_itr++){
      ZeeCategory cat;
      std::cout<<"RooSmearer::GetNumberVariables() says "<<GetNumberVariables()<<std::endl;
      cat.SetNumberVariables(GetNumberVariables());//ZFitter.cpp fixes NumberVariables at runtime. The  information is propagated to RooSmearer, which in turn propagates it to SmearingImporter
      cat.categoryIndex1 = region_ele1_itr - importer._regionList.begin();
      cat.categoryIndex2 = region_ele2_itr - importer._regionList.begin();
      cat.data_events = &(data_events_cache[index]);
      cat.mc_events = &(mc_events_cache[index]);
      
      cat.categoryName1 += *region_ele1_itr;
      cat.categoryName2 += *region_ele2_itr;

      cat.pars1.add(_params_vec[cat.categoryIndex1]);
      cat.pars2.add(_params_vec[cat.categoryIndex2]);

      cat.scale1=1;
      cat.alpha1=0;
      cat.constant1=0;
      cat.scale2=1;
      cat.alpha2=0;
      cat.constant2=0;

      cat.nBins=nBins_;
      cat.targetVariable_min=targetVariable_min_;
      cat.targetVariable_max=targetVariable_max_;
      
      // allocating the histograms
      TString histoName;
      histoName= cat.categoryName1+cat.categoryName2+"_";
      histoName.ReplaceAll("-","_");

      std::cout<<"Setting histo names in RooSmearer::InitCategories"<<std::endl;
      std::cout<<cat.categoryName1<<" "<<cat.categoryName2<<std::endl;
      for(int var=0;var<cat.GetNumberVariables();var++){//Each category can have more than 1 targetVariable
	TString varName="";
	varName.Form("%d",var);
	cat.hist_mc[var]=new TH1F(histoName+varName+"_MC", histoName+varName+"_MC", cat.nBins[var], cat.targetVariable_min[var], cat.targetVariable_max[var]);
	(cat.hist_mc[var])->Sumw2();
	cat.smearHist_mc[var]=new TH1F(histoName+varName+"_smearMC", histoName+varName+"_smearMC", cat.nBins[var], cat.targetVariable_min[var], cat.targetVariable_max[var]);
	(cat.smearHist_mc[var])->Sumw2();
	cat.hist_data[var]=new TH1F(histoName+varName+"_data", histoName+varName+"_data", cat.nBins[var], cat.targetVariable_min[var], cat.targetVariable_max[var]);
	(cat.hist_data[var])->Sumw2();
	cat.smearHist_data[var]=new TH1F(histoName+varName+"_smeardata", histoName+varName+"_smeardata", cat.nBins[var], cat.targetVariable_min[var], cat.targetVariable_max[var]);
	(cat.smearHist_data[var])->Sumw2();
      }//loop over targetVariables for category

      std::cout<<"SetHisto in RooSmearer"<<std::endl;
      SetHisto(*(cat.mc_events), cat.hist_mc);
      SetHisto(*(cat.data_events), cat.hist_data);
      if(_isDataSmeared){//you don't need bool mcToy anymore
	if(cat.data_events->size()!=0){
	  std::cout<<"[INFO toyMC]: In the runToy scenario you want to apply scale/sigma to your data"<<std::endl;
	  //cat.hist_data[0]->SaveAs("tmp/before"+cat.categoryName1+cat.categoryName2+".root");//
	  SetSmearedHisto(*(cat.data_events),
			  cat.pars1, cat.pars2, 
			  cat.categoryName1, cat.categoryName2, 1,
			  cat.smearHist_data);
	}
      }
      cat.smearHist_data[0]->SaveAs("tmp/after"+cat.categoryName1+cat.categoryName2+".root");

      std::cout<<"[src/RooSmearer.cc] Check if the category must be deactivated"<<std::endl;
      //------------------------------ deactivating category in case a targetVariable has few events/bad shape (a shoulder)
      cat.active=true;
      //LOOP DEACTIVATING CATEGORY
      //unsigned int deactiveMinEvents = 1000; 
      //unsigned int deactiveMinEvents = 1;
      //r(int var=0;var<cat.GetNumberVariables();var++){
      //if(cat.categoryIndex1 != cat.categoryIndex2) // not diagonal category
      //  //deactiveMinEvents=1500; //_deactive_minEventsOffDiag;
      //if(cat.hist_mc[var]->Integral() < deactiveMinEvents){
      //  std::cout<<"[INFO] Not enough MC!: "<<cat.hist_mc[var]->Integral()<<std::endl;
      //  std::cout << "[INFO] Category: " << ZeeCategories.size() 
      //	    << ": " << cat.categoryName1 << "\t" << cat.categoryName2
      //	    << " has been deactivated (nEvents < "; 
      //  if(cat.categoryIndex1 != cat.categoryIndex2) std::cout << _deactive_minEventsOffDiag << ")";
      //  else std::cout << _deactive_minEventsDiag << ")";
      //  std::cout << " nEvents mc targetVariable: " << cat.hist_mc[var]->Integral()
      //	    << std::endl;
      //  cat.active=false;
      //}     

      //shoulder check
//	float max=(cat.hist_mc[var])->GetMaximum();
//	float left=(cat.hist_mc[var])->GetBinContent(1);
//	float right=(cat.hist_mc[var])->GetBinContent((cat.hist_mc[var])->GetNbinsX());
//	//if(targetVariable_=="invMass"){
//	//deactivate for high shoulder (irregular shape)
//	if((right - left)/max > 0.2 || (left - right)/max > 0.4){
//	  cat.active=false;
//	  std::cout << "[INFO] Category: " << ZeeCategories.size() 
//		    << ": " << cat.categoryName1 << "\t" << cat.categoryName2
//		    << " has been deactivated for high shoulder (in the MC distribution): " << right << " " << left << " " << max 
//		    << std::endl;
//	}
//      }//end loop deactivate cat

      if (_autoBin){
	AutoNBins(cat);
      }
      cat.nSmearToy=_nSmearToy;
      if(smearscan || (cat.active && _autoNsmear)) AutoNSmear(cat);
      
      cat.nLLtoy=1;
      cat.nll=0;
      std::cout << "[INFO] Category: " << ZeeCategories.size()<<std::endl 
		<< "\t" << cat.categoryName1 << "\t" << cat.categoryName2<<std::endl
		<< "\t" << cat.categoryIndex1 << "\t" << cat.categoryIndex2<<std::endl
		<< "\t" << cat.nSmearToy <<std::endl
		<< "\t" << cat.nLLtoy<<std::endl
		<< "\t" << cat.mc_events->size()<<std::endl 
		<< "\t" << cat.data_events->size()<<std::endl
		<< "\t" << cat.active<<std::endl;
      for(int var=0;(unsigned)var <cat.GetNumberVariables();var++){
	  std::cout<< "\t" << cat.nBins[var] <<std::endl
		   << "\t" << (cat.hist_mc[var])->Integral()<<std::endl
		   << "\t" << (cat.hist_data[var])->Integral()<<std::endl
	    //<< "\t" << right[var]
	    //<< "\t" << left[var]
	    //<< "\t" << max[var]
		   << std::endl;
      }
      cat.hist_data[0]->SaveAs("tmp/after_justBeforePush_back"+cat.categoryName1+cat.categoryName2+".root");
      ZeeCategories.push_back(cat);      
      index++;	 
    }
  }
}


std::vector<TH1F *> RooSmearer::GetSmearedHisto(ZeeCategory& category, bool isMC,
				  bool smearEnergy, bool forceNew, bool multiSmearToy){
  //std::cout<<"[STATUS] Calling RooSmearer::GetSmearedHisto"<<std::endl;
  std::vector<TH1F *> histo_vector;//final vector to be returned
  std::vector<TH1F *> h;//temp vector
  
  for(int var=0;var<category.GetNumberVariables();var++){
    histo_vector.push_back(NULL);
    h.push_back(NULL);
  }
  
  if(isMC){
    h = (smearEnergy) ? category.smearHist_mc : category.hist_mc;
  }else{
    //h = (smearEnergy) ? category.smearHist_data : category.hist_data;//This is the original
    if(smearEnergy){
      h=category.smearHist_data;
    }else{
      h=category.hist_data;
      //std::cout<<"for data h is category.hist_data"<<std::endl;
    }
  }
  
  for(int var=0;var<category.GetNumberVariables();var++){
    if(h[var]==NULL){
      std::cerr << "[ERROR] histogram pointer in RooSmearer::GetSmearedHisto is NULL!!! for var " <<var<<std::endl;
      exit(1);
    }
  }
  
  //Taking the right histograms
  if(smearEnergy==false /*&& h->GetEntries() != 0*/){return histo_vector=h;}
  if(smearEnergy==true && isMC==false /*&& h->GetEntries() != 0*/){return histo_vector=h;}

  //You want to smear to MC, not the data
  zee_events_t *cache = (isMC) ? category.mc_events : category.data_events;
  if(cache->size()==0){return histo_vector=h;}
  
  if(smearEnergy){
    for(int var=0;var<category.GetNumberVariables();var++){
      h[var]->Reset();
    }
    
    if(isMC && multiSmearToy){
      SetSmearedHisto(*cache,
		      category.pars1, category.pars2, 
		      category.categoryName1, category.categoryName2, category.nSmearToy,
		      h);
    } else {
      SetSmearedHisto(*cache,
		      category.pars1, category.pars2, 
		      category.categoryName1, category.categoryName2, 1,
		      h);
    }
    
    if(isMC){//! normalize mc histo, to be a pdf
      for(int var=0;var<category.GetNumberVariables();var++){
	//std::cout<<"mc normalization in RooSmearer::SetSmearedHisto is "<<h[var]->Integral()<<std::endl;
	double int_mc=h[var]->Integral();
	h[var]->Scale(1./h[var]->Integral());//Integral avoids underflows and overflows (I checked it by hand)
	h[var]->SetBinContent(0,int_mc);//the underflow contains the original MC normalization
	//h[var]->Scale(1./h[var]->Integral(1,h[var]->GetNbinsX()));
      }
    }
    //(*h)->Smooth();//I removed the smoothing
  } else{
    for(int var=0;var<category.GetNumberVariables();var++){
      h[var]->Reset(); // unuseful: histogram should be empty
    }
    SetHisto(*cache, h);
  }
  
  histo_vector=h;
  return histo_vector;
}

void RooSmearer::SetHisto(const zee_events_t& cache, std::vector<TH1F *> hist) const{
  //std::cout<<"[STATUS] Calling RooSmearer::SetHisto"<<std::endl;
  for(int var=0;var<hist.size();var++){
    hist[var]->Reset();
  }

  for(zee_events_t::const_iterator event_itr = cache.begin(); 
      event_itr!= cache.end();
      event_itr++){
    for(int var=0;var<hist.size();var++){
      hist[var]->Fill(event_itr->targetVariable[var],event_itr->weight);
      //if(event_itr->targetVariable[1]<40){
      //std::cout<<"hei, pt 2 is "<<event_itr->pt2<<std::endl;
      //std::cout<<"hei, pt 1 is "<<event_itr->pt1<<std::endl;
      //std::cout<<"hei, energy 2 is "<<event_itr->energy_ele2<<std::endl;
      //std::cout<<"hei, energy 1 is "<<event_itr->energy_ele1<<std::endl;
      //
      //}
    }
  }

  return;
}

void RooSmearer::SetSmearedHisto(const zee_events_t& cache, 
				 RooArgSet pars1, RooArgSet pars2, 
				 TString categoryName1, TString categoryName2, unsigned int nSmearToy,
				 std::vector<TH1F *> hist) const{
  //std::cout<<"Inside RooSmearer::SetSmearedHisto"<<std::endl;
  //std::cout<<"Target Variable attribute in RooSmearer is "<<targetVariable_<<std::endl;
  // retrieve values from params for the category
  float scale1 = pars1.getRealValue("scale_"+categoryName1, 0., kTRUE);
  float constant1 = pars1.getRealValue("constTerm_"+categoryName1, 0., kTRUE);
  float alpha1 = pars1.getRealValue("alpha_"+categoryName1, 0., kTRUE);

  float scale2 = pars2.getRealValue("scale_"+categoryName2, 0., kTRUE);
  float constant2 = pars2.getRealValue("constTerm_"+categoryName2, 0., kTRUE);
  float alpha2 = pars2.getRealValue("alpha_"+categoryName2, 0., kTRUE);

#ifdef DEBUG
  pars1.writeToStream(std::cout, kFALSE);
  std::cout << "--" << std::endl;
  pars2.writeToStream(std::cout, kFALSE);
#endif

  double smearEne1[1000], smearEne2[1000]; // _nSmearToy<100
  if(cache.begin()->smearings_ele1==NULL){
    std::cerr<< "[ERROR] No smearings" << std::endl;
    exit(1);
  }
  int i=0;
  for(zee_events_t::const_iterator event_itr = cache.begin(); 
	  event_itr!= cache.end();
	  event_itr++){

    smearedEnergy(smearEne1, nSmearToy, event_itr->energy_ele1, scale1, alpha1, constant1, event_itr->smearings_ele1);
    smearedEnergy(smearEne2, nSmearToy, event_itr->energy_ele2, scale2, alpha2, constant2, event_itr->smearings_ele2);
    //std::cout<<"In SetSmearedHisto hist.size is "<<hist.size()<<std::endl;
    for(unsigned int iSmearToy=0; iSmearToy < nSmearToy; iSmearToy++){
      for(int var=0;var<hist.size();var++){//hist is a vector<TH1F *>
	if(targetVariable_[var]=="invMass"){
	  hist[var]->Fill(event_itr->targetVariable[var] * sqrt(smearEne1[iSmearToy] * smearEne2[iSmearToy]),event_itr->weight);
	}else if(targetVariable_[var]=="ptRatio"){
	  if(!(event_itr->isSwapped)){
	    hist[var]->Fill(event_itr->targetVariable[var] * (smearEne1[iSmearToy]/smearEne2[iSmearToy]), event_itr->weight);
	  }else{
	    hist[var]->Fill(event_itr->targetVariable[var] * (smearEne2[iSmearToy]/smearEne1[iSmearToy]), event_itr->weight);
	  }
	}else if(targetVariable_[var]=="ptSum"){
	  if(!(event_itr->isSwapped)){
	    hist[var]->Fill((smearEne1[iSmearToy]*event_itr->pt1 + smearEne2[iSmearToy]*event_itr->pt2), event_itr->weight);
	  }else{
	    hist[var]->Fill((smearEne1[iSmearToy]*event_itr->pt2 + smearEne2[iSmearToy]*event_itr->pt1), event_itr->weight);
	  }
	}else if(targetVariable_[var]=="pt2Sum"){
	  if(!(event_itr->isSwapped)){
	    //std::cout<<"target "<<event_itr->targetVariable[var]<<std::endl;
	    //std::cout<<"pt1 "<<event_itr->pt1<<std::endl;
	    //std::cout<<"pt2 "<<event_itr->pt2<<std::endl;
	    //std::cout<<"smear1 "<<smearEne1[iSmearToy]<<std::endl;
	    //std::cout<<"smear2 "<<smearEne2[iSmearToy]<<std::endl;
	    //double smeared_variable=sqrt(pow(smearEne1[iSmearToy]*event_itr->pt1,2) + pow(smearEne2[iSmearToy]*event_itr->pt2,2));
	    //std::cout<<"smeared variable "<<smeared_variable<<std::endl;
	    hist[var]->Fill(sqrt(pow(smearEne1[iSmearToy]*event_itr->pt1,2) + pow(smearEne2[iSmearToy]*event_itr->pt2,2)), event_itr->weight);
	  }else{
	    hist[var]->Fill(sqrt(pow(smearEne1[iSmearToy]*event_itr->pt2,2) + pow(smearEne2[iSmearToy]*event_itr->pt1,2)), event_itr->weight);
	  }
	}	        
      }//loop over vector<TH1> (more than 1 variable per category)
    }//loop over toys
  }//loop over cache events

  for(int var=0;var<hist.size();var++){
    hist[var]->Scale(1./nSmearToy);
  }
  return;
}



double RooSmearer::smearedEnergy(double *smear, unsigned int nGen, float ene,float scale,float alpha,float constant, const float *fixedSmearings) const
{
  // sigmaMB = sigma Material Budget
  // if I want to take into account the non perfet simulation of the
  // material budget in front of ECAL, I have to add in quadrature an
  // additional smearing term sigmaMB

  // float sigmaMB = sigmaMB(ene, eta, phi);
  // sigmaMB non dovrebbe essere una smooth function, ma tenere in
  // conto differenze consistenti in piccole regioni.
  //  float sigmaMB = MB * fabs{etaEle} * fabs{etaEle}
  //   float sigma = sqrt(alpha/sqrt(ene)*alpha/sqrt(ene) + constant * constant + sigmaMB*sigmaMB);

  float sigma = sqrt(alpha*alpha/ene + constant * constant );
  if(sigma==0){
    for(unsigned int i=0; i<nGen; i++){
      smear[i] = scale;
    }
  } else {
#ifdef FIXEDSMEARINGS
  for(unsigned int i=0; i < NSMEARTOYLIM && i<nGen; i++){
    smear[i] = (double) (fixedSmearings[i]*sigma)+(scale);
  }
  for(unsigned int i=NSMEARTOYLIM; i < nGen; i++){
    smear[i] = rgen_->Gaus(scale,sigma);
  }
#else
  for(unsigned int i=0; i < nGen; i++){
	smear[i] = rgen_->Gaus(scale,sigma);
  }
#endif
  }
  return smear[0];
}



double RooSmearer::getLogLikelihood(std::vector<TH1F*> data, std::vector<TH1F*> prob) const
{//prob is the normalized mc
  //it returns the loglikelihood PER CATEGORY (then you should minimize the sum of the NLLs)
  for(int var=0;var<data.size();var++){
    if (!data[var] || !prob[var])
      {
	std::cout << "ERROR: at least one empty histograms given in the vector" << std::endl;
	return -9999999.;
      }
  }

  for(int var=0;var<data.size();var++){
    if (data[var]->GetNbinsX() != prob[var]->GetNbinsX())
      {
	std::cout << "ERROR: data and probabilities for the histogram # "<<var<<" are binned differently" << std::endl;
	return -9999999.;
      }
  }

  double logL=0.;
  //int flag_empty=0;
  //Not using underflows and overflows at the moment
  for(int var=0;var<data.size();var++){
    for (int ibin=1;ibin<data[var]->GetNbinsX()+1;++ibin)
      {
	
	if(prob[var]->GetBinContent(ibin)==0 && data[var]->GetBinContent(ibin)!=0){
	  //PENALTY TERM
	  //flag_empty=1;
	  logL+=data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(1./(prob[var]->GetBinContent(0)));//underflows contains the MC normalization (initial)
	  //std::cout<<"probabilita' era zero ****"<<std::endl;
	  //std::cout<<"normalization "<<prob[var]->GetBinContent(0)<<std::endl;
	  //std::cout<<"1 over is "<<(1./prob[var]->GetBinContent(0))<<std::endl;
	  //std::cout<<"contribution is "<<data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(1./(prob[var]->GetBinContent(0)))<<std::endl;
	}else if(prob[var]->GetBinContent(ibin)>0){//less than zero if possible for mcGenWeight<0, but then ln(p) is a random rumber
	  //double weight= 	(data[var]->GetBinContent(ibin)>0) ? data[var]->GetBinError(ibin) * data[var]->GetBinError(ibin) / data[var]->GetBinContent(ibin) : 1;
	  double weight= 1.;
	//	std::cout << "weight = " << weight << std::endl;
#ifdef POISSON_LIKELIHOOD
	  //Poisson likelihood 	  
	  logL+=  weight * (data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(prob[var]->GetBinContent(ibin)) - prob[var]->GetBinContent(ibin));
#else
	// Binomial likelihood
	  //std::cout<<"** bin "<<ibin<<" "<<weight * (data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(prob[var]->GetBinContent(ibin)) + (data[var]->Integral()-data[var]->GetBinContent(ibin))*ROOT::Math::Util::EvalLog(1-prob[var]->GetBinContent(ibin)))<<std::endl;
//logL+= weight * (data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(prob[var]->GetBinContent(ibin)) + (data[var]->Integral()-data[var]->GetBinContent(ibin))*ROOT::Math::Util::EvalLog(1-prob[var]->GetBinContent(ibin)));

//Multinomial likelihood
	  //std::cout<<"** bin "<<ibin<<" n_i: "<<data[var]->GetBinContent(ibin)<<" p_i: "<<prob[var]->GetBinContent(ibin)<<" ln(p_i): "<<ROOT::Math::Util::EvalLog(prob[var]->GetBinContent(ibin))<<" lnL: "<<data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(prob[var]->GetBinContent(ibin))<<std::endl;
	  logL+=data[var]->GetBinContent(ibin)* ROOT::Math::Util::EvalLog(prob[var]->GetBinContent(ibin));
#endif
	} 
      }
  }// loop over var for each varible of the targetVariable
  //std::cout<<"In RooSmearer::getLogLikelihood, the negative loglikelihood of the specific category is "<<-logL<<std::endl;

  //if(flag_empty==1){
    //logL*=1.5;//penalty
    //std::cout<<"with penalty, -logL is "<<-logL<<std::endl;
  //}
  return logL;
}



double RooSmearer::getCompatibility(bool forceUpdate) const
{
  //std::cout<<"[STATUS] In RooSmearer::getCompatibility returns the total NLL (negative log likelihood): sum of the nll per category"<<std::endl;
  RooSmearer* myClass=(RooSmearer *) this;

  std::vector<TH1F *> mcHistos, dataHistos;

  double compatibility = 0.;
  bool updated=false;
  for(std::vector<ZeeCategory>::iterator cat_itr = myClass->ZeeCategories.begin();
      cat_itr != myClass->ZeeCategories.end();
      cat_itr++){
    if(!cat_itr->active) continue;
    if(forceUpdate||isCategoryChanged(*cat_itr,true)){ // && withSmearToy){
      updated=true;
      #ifdef DEBUG
      std::cout << "[DEBUG] " << cat_itr->categoryName1 << " - " << cat_itr->categoryName2 << "\t isupdated" << std::endl;
      #endif
      myClass->UpdateCategoryNLL(*cat_itr, cat_itr->nLLtoy); //the new nll has been updated for the category
    }
    compatibility+=cat_itr->nll;
    myClass->lastNLLrms+= (cat_itr->nllRMS*cat_itr->nllRMS);
    //    std::cout << std::fixed << std::setprecision(20) << "[DEBUG] Compatibility: " << compatibility << "\t" << compatibility-lastNLL << std::endl;
  }
  //if(withSmearToy) std::cout << "[DEBUG] Compatibility2: " << compatibility << "\t" << compatibility - nllMin << std::endl;
  if(dataset!=NULL && updated){
    myClass->nllVar.setVal(compatibility);
    myClass->dataset->add(RooArgSet(_paramSet,nllVar));
  }

  myClass->lastNLL=compatibility;
  myClass->lastNLLrms=sqrt(lastNLLrms);
  
  if(nllBase==0) myClass->nllBase=compatibility;
  else {
    if(nllMin> compatibility){
      //  std::cout << "[DEBUG] update nllMin: " << nllMin << " -> " << compatibility << std::endl;
      myClass->nllMin=compatibility;
    }
  }
  #ifdef DEBUG
  std::cout<<"updated chi2 (from RooSmearer::getCompatibility) is "<<compatibility<<std::endl;
  #endif
  return compatibility;
}



Double_t RooSmearer::evaluate() const
{
  //evaluate gives the chi2: i.e. the minimization parameter. The Chi2 is the sum of the NLL of the categories
  //update last result
  double comp_mean = getCompatibility();//comp_mean IS the chi2
  //std::cout<<"Chi2 is "<<comp_mean<<std::endl;
  RooSmearer* myClass=(RooSmearer *) this;
  double weight=(nllBase*2-comp_mean);
  if(weight<0) weight=1;
  myClass->_markov.AddFast(myClass->_paramSet, comp_mean, weight);
  return comp_mean;
}




int RooSmearer::Trag_eq(int row, int col, int N) const
{
  if (row <= col)
	return row*N - (row-1)*((row-1) + 1)/2 + col - row;
  else 
	return col*N - (col-1)*((col-1) + 1)/2 + row - col;
}



/*NOT USED
void RooSmearer::SetAutoBin(ZeeCategory& category,double min, double max){
  std::cout<<"Calling RooSmearer::SetAutoBin"<<std::endl;
  category.hist_mc[0]->SetBins(1000,min,max);
  category.hist_mc[0]->Reset();
  TH1F* wide = GetSmearedHisto(category,true, false);
  Int_t nq = 100;
  double xq[100],yq[100];
  for (Int_t i=0;i<nq;i++) 
    xq[i] = (double)(i+1)/nq;
  wide->GetQuantiles(nq,yq,xq);
  

  // center the mean of the distribution
  float middle=0.5*(category.targetVariable_max+category.targetVariable_min);
  
  if(fabs(middle - wide->GetMean()) > 5){
    int shift = (int)(wide->GetMean() - middle);
    std::cout << "[INFO] Shifting targetVariable range of " << shift << std::endl;
    std::cout << "       targetVariable_min = " << category.targetVariable_min+shift << std::endl;
    std::cout << "       targetVariable_max = " << category.targetVariable_max+shift << std::endl;
    category.targetVariable_min+=shift;
    category.targetVariable_max+=shift;
  }
  ResetBinning(category);
}
*/

void RooSmearer::ResetBinning(ZeeCategory& category){
  std::cout<<"Calling RooSmearer::ResetBinning"<<std::endl;
  for(int var=0;var<category.GetNumberVariables();var++){
    category.hist_mc[var]->SetBins(category.nBins[var], category.targetVariable_min[var], category.targetVariable_max[var]);
    category.hist_mc[var]->Reset();
    category.smearHist_mc[var]->SetBins(category.nBins[var], category.targetVariable_min[var], category.targetVariable_max[var]);
    category.smearHist_mc[var]->Reset();
    
    category.hist_data[var]->SetBins(category.nBins[var], category.targetVariable_min[var], category.targetVariable_max[var]);
    category.hist_data[var]->Reset();
    if(category.smearHist_data[var]!=NULL){
      category.smearHist_data[var]->SetBins(category.nBins[var], category.targetVariable_min[var], category.targetVariable_max[var]);
      category.smearHist_data[var]->Reset();
    }
  }
  return;
}


void RooSmearer::AutoNBins(ZeeCategory& category){
  std::cout<<"Calling RooSmearer::AutoNBins"<<std::endl;
  if(!category.active) return;

  for(int var=0; var<category.GetNumberVariables();var++){
    if(category.hist_mc[var]->Integral()>30000){
      std::cout << "[INFO] Category: " << category.categoryName1 << " " << category.categoryName2 << "\t binning increased to: "  <<     category.nBins[var]*2 << std::endl;
      category.nBins[var]*=2;
      ResetBinning(category);
    }else   if(category.hist_mc[var]->Integral()<5000){
      std::cout << "[INFO] Category: " << category.categoryName1 << " " << category.categoryName2 << "\t binning reduced to: "  <<     category.nBins[var]/2 << std::endl;
      category.nBins[var]/=2;
      ResetBinning(category);
    }
  }
  
  SetHisto(*(category.mc_events),   category.hist_mc);
  SetHisto(*(category.data_events), category.hist_data);
  if(_isDataSmeared){
    std::cout<<"[INFO toyMC]: Rebinning SmearHist Data"<<std::endl;
    SetSmearedHisto(*(category.data_events),
		    category.pars1, category.pars2, 
		    category.categoryName1, category.categoryName2, 1,
		    category.smearHist_data);
  }
  

  return;

    /* MA VIENE CHIUSO PRIMA A CHE SERVE??
    //------------------------------ rescale mc histograms 
  
    //  (isCategoryChanged(category,true)); // && withSmearToy){
    category.scale1=1;
    TH1F *data = GetSmearedHisto(category, false, _isDataSmeared, true, false); ///-----> not need to repeate!
    TH1F *mc = GetSmearedHisto(category, true, true,true); // regenerate the histogram: forceNew
    double chi2old=getLogLikelihood(data, mc);
    double min=chi2old; int nBin_min=0;

    std::cout << "[AUTOBIN] Starting scale1 sensitivity check: " << category.scale1 << "\t" <<std::fixed << std::setprecision(10) << chi2old << std::endl;
    double valueMin=category.scale1;

    for(int n=0; n<10; n++){
      double value=category.scale1-0.01+0.002*n;
      category.pars1.setRealValue("scale_"+category.categoryName1, value , kTRUE);

      ResetBinning(category);
      data = GetSmearedHisto(category, false, _isDataSmeared, true, false); ///-----> not need to repeate!
      mc = GetSmearedHisto(category, true, true,true); // regenerate the histogram: forceNew
      double chi2=getLogLikelihood(data, mc);
      if(chi2 < min){
	min=chi2;
	valueMin=value;
      }
      std::cout << "[BINNING] scale1: " << value << "\t" << chi2-chi2old << "\t" << chi2-min << std::endl;
    }
    if(category.scale1!=valueMin){
      std::cout << "[WARNING] scale1 not closing: " << category.scale1 << "\t" << valueMin << std::endl;
    }
    return;

    category.nBins=nBin_min;
    ResetBinning(category);
    return;
    */

}

void RooSmearer::AutoNSmear(ZeeCategory& category){
  //??
  
  if(!category.active) return;
  //AutoNBins(category);
  return;
  //  category.nSmearToy=std::max(NSMEARTOYLIM*2,200);
  //Long64_t catSize=category.mc_events->size();
//   if(catSize<3000){ 
//     category.nBins=40;
//     ResetBinning(category);
//   } else if(catSize>40000){
//     category.nBins=160;
//     ResetBinning(category);
//   }
  
//   category.nSmearToy=(int)(300000./catSize); 
//   if(category.nSmearToy <7) category.nSmearToy = 7; // fix the min to 3
//   else  if( category.nSmearToy > 40) category.nSmearToy = 40; // fix the max to 20
  //  category.nSmearToy=NSMEARTOYLIM;
  if(!smearscan) return;

  //------------------------------ rescale mc histograms 
  double stdDev=10, stdDevLim=0.3;
  double min=10; //int nBin_min=0;
  int n=0;
  unsigned int nSmearToyLim=1000;
  category.nLLtoy=1;
  for(int iBin=160; iBin>120; iBin/=2){
    //category.nBins= iBin;
    //ResetBinning(category);
    //TH1F *data = GetSmearedHisto(category, false, _isDataSmeared, false, false); ///-----> not need to repeate!
    GetSmearedHisto(category, false, _isDataSmeared, false, false); ///-----> not need to repeate!

    for(category.nLLtoy=1; category.nLLtoy < 2; category.nLLtoy+=2){
      for(; category.nSmearToy <= nSmearToyLim && stdDev> stdDevLim; category.nSmearToy*=2){
	double  sum=0, sum2=0;
	TStopwatch cl;
	cl.Start();
	for(n=0; n<50; n++){
	  UpdateCategoryNLL(category, category.nLLtoy); //the new nll has been updated for the category
	  
	  sum+=category.nll;
	  sum2+=category.nll * category.nll;
	}
	cl.Stop();
	sum/=n;
	sum2/=n;
	stdDev= sqrt(sum2 - sum*sum);
	if(stdDev/sum < min) {
	  min=stdDev/sum; 
	  //nBin_min=category.nBins;
	}

	
	std::cout << "[SMEARSCAN] " 
		  << "\t" << category.categoryIndex1
		  << " " << category.categoryIndex2
		  << "\t" << stdDev 
		  << "\t" << sum
		  << " " << stdDev/sum 
		  << "\t" << category.nSmearToy 
		  << " " << category.nLLtoy
	  //<< " " << category.nBins[var] 
		  << "\tmcEvents=" << category.mc_events->size() 
		  << "\tdataEvents=" << category.data_events->size() 
		  << "\t" << cl.CpuTime()/50. << std::endl; 
	//<< std::endl;
      }
      //if(category.nSmearToy > nSmearToyLim) category.nSmearToy/=2;
      category.nSmearToy/=2;
    }
  }
  
  //  exit(0);
  return;
#ifdef DD
  //============================================================

  AutoNBins(category);
  
  Long64_t catSize=category.mc_events->size();
  category.nSmearToy=(int)(std::max(1.,1000000./catSize));
  std::cout 
    << "\t" << category.categoryName1 << " " << category.categoryName2 
    << "\tnSmearToy=" << category.nSmearToy 
    << std::endl;
  return; 
#ifndef DUMPHIST
  //  AutoNBins(category);
#endif
  if(!category.active) return;
  TH1F *data = GetSmearedHisto(category, false, _isDataSmeared, false, false); ///-----> not need to repeate!


  //------------------------------ rescale mc histograms 
  double stdDev=0, mean=0;
  int n=0;


  float nSmearToyLim=std::max(1.,100000./catSize);
  for(category.nLLtoy=1; category.nLLtoy <= 1; category.nLLtoy+=2){
      for(category.nSmearToy=1; category.nSmearToy <= nSmearToyLim; category.nSmearToy+=5){

    double  sum=0, sum2=0;
    for(n=0; n<500; n++){
      double comp=0, comp2=0.;
      for(unsigned int itoy=0; itoy < category.nLLtoy; itoy++){
	TH1F *mc = GetSmearedHisto(category, true, true,true); // regenerate the histogram: forceNew
#ifdef DUMPHIST
	TString histName="tmp/hist3/";
	histName+=category.mc_events->size();
	histName+="-"; histName+=category.nSmearToy;
	histName+="-"; histName+=category.nLLtoy;
	histName+="-"; histName+=category.nBins;
	histName+="-"; histName+=n;
	histName+="-"; histName+=importer._constTermToy;
	histName+=".root";
	mc->SaveAs(histName);
#endif   
	double c = getLogLikelihood(data, mc);
	comp+=c;
	comp2+=c*c;
      }
      comp/=category.nLLtoy;
      comp2/=category.nLLtoy;
      sum+=comp;
      sum2+=comp*comp;
    }
    
    mean=sum/n;
    stdDev= sqrt(sum2/n - mean*mean);
    std::cout << "[DEBUG] n = " << n 
	      << "\t" << category.categoryIndex1
	      << "\t" << category.categoryIndex2
	      << "\t" << stdDev 
	      << "\t" << mean 
	      << "\t" << stdDev/mean 
	      << "\t" << category.nSmearToy 
	      << "\t" << category.nLLtoy
	      << "\t" << category.nBins 
	      << "\tmcEvents=" << category.mc_events->size() 
	      << "\tdataEvents=" << category.data_events->size() 
	      << std::endl;
      }
  }


  exit(0);

  do{
    category.nSmearToy+=5;
    double  sum=0, sum2=0;
    for(n=0; n<500; n++){
      TH1F *mc = GetSmearedHisto(category, true, true, true); // forceNew
      TH1F *h = (TH1F *) mc->Clone(TString(mc->GetName())+"_c");
      h->Scale(1./h->Integral()); // now the histogram is a pdf 
      double compatibility = fabs(getLogLikelihood(data, h ));
      sum+=compatibility;
      sum2+=compatibility *compatibility;
      stdDev=sum2/(n+1) - (sum/(n+1))* (sum/(n+1));
      //std::cout << "[DEBUG] n = " << n << "\t" << compatibility << "\t" << stdDev << "\t" << sum2/(n+1) << "\t" << sum/(n+1) << std::endl;
    delete h;
    }
    mean=sum/n;
    stdDev= sqrt(sum2/n - mean*mean);
    std::cout << "[DEBUG] n = " << n 
	      << "\t" << stdDev 
	      << "\t" << mean 
	      << "\t" << stdDev/mean 
	      << "\t" << category.nSmearToy 
	      << "\t" << category.nBins 
	      << "\tmcEvents=" << category.mc_events->size() 
	      << "\tdataEvents=" << category.data_events->size() 
	      << std::endl;
  } while (stdDev/mean>0.0005 && category.nSmearToy<2);
  std::cout << "[INFO] Category: " 
	    << category.categoryIndex1 << " " << category.categoryIndex2 
    //<< "\t" << category.categoryName1 << " " << category.categoryName2 
	    << "\tstdDev=" << stdDev << "\tnSmearToy="<< category.nSmearToy 
	    << "\tmcEvents=" << category.mc_events->size() 
	    << "\tdataEvents=" << category.data_events->size() 
	    << std::endl;

  //TH1F *mc = GetSmearedHisto(category, true, true); // forceNew
  //  data->SaveAs("data.root");
  //mc->SaveAs("mc.root");
#endif
  exit(0);
  return;
}



bool RooSmearer::isCategoryChanged(ZeeCategory& category, bool updateVar) const{
#ifdef DEBUG
  std::cout << "[DEBUG] Checking if category changed: " << category.categoryName1 << "\t" << category.categoryName2 << std::endl;
#endif
    bool changed=false;

    RooArgList argList1(category.pars1);
    RooArgList argList2(category.pars2);

    // checking if one of the variables has changed
    TIterator *it = argList1.createIterator();
    for(RooRealVar *v = (RooRealVar *) it->Next(); v!=NULL; v = (RooRealVar*) it->Next()){
      TString varName=v->GetName();
      double varValue=v->getVal();
      if(varName.Contains("scale") && varValue!= category.scale1){
	changed=true;
#ifdef DEBUG
	std::cout << "scale changed for: " << category.categoryName1 << "\t" << category.scale1 << "\t" << varValue << "\t" << changed << "\t" << updateVar << std::endl;
#endif
	if(updateVar) category.scale1=varValue;
      }
      if(varName.Contains("alpha")     && varValue!= category.alpha1){
	changed=true;
#ifdef DEBUG
	std::cout << "alpha changed for: " << category.categoryName1 << "\t" << category.alpha1 << "\t" << varValue << "\t" << changed << std::endl;
#endif
	if(updateVar) category.alpha1=varValue;
      }
      if(varName.Contains("constTerm") && varValue!= category.constant1){
	changed=true;
#ifdef DEBUG
	std::cout << "alpha changed for: " << category.categoryName1 << "\t" << category.constant1 << "\t" << varValue << "\t" << changed << std::endl;
#endif
	if(updateVar) category.constant1 = varValue;
      }
    }

    delete it;
    it = argList2.createIterator();
    for(RooRealVar *v = (RooRealVar *) it->Next(); v!=NULL; v = (RooRealVar*) it->Next()){
      TString varName=v->GetName();
      double varValue=v->getVal();

      if(varName.Contains("scale")     && varValue!= category.scale2){
	changed=true;
	if(updateVar) category.scale2=varValue;
      }
      if(varName.Contains("alpha")     && varValue!= category.alpha2){
	changed=true;
	if(updateVar) category.alpha2=varValue;
      }
      if(varName.Contains("constTerm") && varValue!= category.constant2){
	changed=true;
	if(updateVar) category.constant2 = varValue;
      }
    }
    
    delete it;
    return changed;
}      


void RooSmearer::SetNSmear(unsigned int n_smear, unsigned int nlltoy){

  for(std::vector<ZeeCategory>::iterator cat_itr = ZeeCategories.begin();
      cat_itr != ZeeCategories.end();
      cat_itr++){
    if(n_smear>0) cat_itr->nSmearToy=n_smear;
    if(nlltoy >0) cat_itr->nLLtoy=nlltoy;
  }
  return;
}


void RooSmearer::Init(TString commonCut, TString eleID, Long64_t nEvents, bool mcToy, bool externToy, TString initFile){
  std::cout<<"[STATUS] In RooSmearer::Init"<<std::endl;
  std::cout<<"[INFO] in RooSmearer::Init nEvents is "<<nEvents<<std::endl;
  std::cout<<"[INfo] commonCut is "<<commonCut<<std::endl;
  std::cout<<"mcToy is "<<mcToy<<std::endl;
  if(mcToy) _isDataSmeared=!externToy; //mcToy;//externToy means that you have a root file with the corrections you want to apply
  std::cout<<"externToy is "<<externToy<<std::endl;
  std::cout<<"_isDataSmeared is "<<_isDataSmeared<<std::endl;
  if(initFile.Sizeof()>1){
    std::cout << "[toyMC INFO] In RooSmearer::Init, truth values for toys initialized to " << std::endl;
    std::cout << "------------------------------ Read init toy MC:" << std::endl;
    _paramSet.readFromFile(initFile);
    _paramSet.writeToStream(std::cout, kFALSE);
  }
  SetCommonCut(commonCut); SetEleID(eleID);
  SetCache(nEvents, mcToy, externToy); 
  InitCategories(mcToy);//the bool mcToy was not used in InitCategories: I restored it
  /// Testing histograms after InitCategories
  //TFile *f = new TFile("tmp/outside_init.root", "recreate");
  TFile *f = new TFile("tmp/outside_init.root", "RECREATE");
  f->Print();
  f->cd();

  std::cout<<"in RooSmearer _isDataSmeared is "<<_isDataSmeared<<std::endl;
  for(std::vector<ZeeCategory>::iterator itr= ZeeCategories.begin();
      itr != ZeeCategories.end();
      itr++){
    //if(!itr->active) continue; 
    std::vector<TH1F *> MC = GetSmearedHisto(*itr, true, false);
    std::vector<TH1F *> smearMC = GetSmearedHisto(*itr, true, true);
    std::vector<TH1F *> data = GetSmearedHisto(*itr, false, _isDataSmeared);
    for(int var=0; var<data.size();var++){
      MC[var]->Write();
      smearMC[var]->Write();
      data[var]->Write();
    }
    f->Write();
    }
  f->Close();
  ///
  TFile *f2 = new TFile("tmp/outside_init2.root", "RECREATE");
  f2->Print();
  f2->cd();

  for(std::vector<ZeeCategory>::iterator itr= ZeeCategories.begin();
      itr != ZeeCategories.end();
      itr++){
    //if(!itr->active) continue; 
    std::vector<TH1F *> data = (*itr).hist_data;
    std::vector<TH1F *> smear_data = (*itr).smearHist_data;
    for(int var=0; var<data.size();var++){
      data[var]->Write();
      smear_data[var]->Write();
    }
    f2->Write();
    }
  f2->Close();
  ///

  TStopwatch cl;
  cl.Start();
  evaluate();
  cl.Stop();
  std::cout << "[INFO] In RooSmearer::Init ==> Time for first eval: "<<std::endl;
  cl.Print();
  if(mcToy && false){//se false non e' settato, i parametri iniziali di fit vengono randomizzati o puoi giocarci come vuoi tu
    RooArgList argList(_paramSet);
    TIterator *it = argList.createIterator();
    for(RooRealVar *var = (RooRealVar *) it->Next(); var!=NULL; var =  (RooRealVar *)it->Next()){
      //var->randomize();
      TString name=var->GetName();
      if(name.Contains("scale_EE-invMass_100_2000-Et_25-trigger-noPF") || name.Contains("scale_EB-invMass_100_2000-Et_25-trigger-noPF")){
	var->setVal(0.97);
      }else if (name.Contains("constTerm")){
	var->setVal(0.02);
      }else {//alpha
	var->setVal(0.);
      }

      std::cout<<"name is "<<var->GetName()<<std::endl;
    }
    std::cout << "------------------------------ Initial value for the fit:" << std::endl;
    _paramSet.writeToStream(std::cout, kFALSE);
  }

  // set initial nll values
  getCompatibility(true);
  return;
}

void RooSmearer::UpdateCategoryNLL(ZeeCategory& cat, unsigned int nLLtoy, bool multiSmearToy){
  //std::cout<<"[INFO] in RooSmearer::UpdateCategoryNLL is called"<<std::endl;

  //TH1F *data = GetSmearedHisto(cat, false, _isDataSmeared,true, false); ///-----> not need to repeate! 1 one smearing! otherwise bin errors are wrongly reduced
  std::vector<TH1F *> data = GetSmearedHisto(cat, false, _isDataSmeared,true, false);
  double comp=0., comp2=0.;
  for(unsigned int itoy=0; itoy < nLLtoy; itoy++){
    //TH1F *mc = GetSmearedHisto(cat, true, true,true,true); // regenerate the histogram: forceNew
    std::vector<TH1F *> mc = GetSmearedHisto(cat, true, true,true,true); // regenerate the histogram: forceNew
    
    double c = getLogLikelihood(data, mc);
    comp+=c;
    comp2+=c*c;
  }
  comp/=nLLtoy;
  comp2/=nLLtoy;

  cat.nll= -comp;
  cat.nllRMS= sqrt(comp2-comp*comp);
  return;
}

void RooSmearer::DumpNLL(void) const{
  std::cout << "[DUMP NLL] " << "Cat1\tCat2\tNLL\tNevt mc\tNevt data\tisActive\tNevt mc\tNevt data" << std::endl;

  for(std::vector<ZeeCategory>::const_iterator cat_itr = ZeeCategories.begin();
      cat_itr != ZeeCategories.end();
      cat_itr++){
    if(!cat_itr->active)
      std::cout << "[DUMP NLL] " << std::setprecision(10) 
		<< cat_itr->categoryIndex1 << " " << cat_itr->categoryIndex2 
		<< "\t" << cat_itr->nll 
		<< "\t" << cat_itr->mc_events->size() << "\t" << cat_itr->data_events->size() 
		<< "\t1" 
		<< "\t" << cat_itr->hist_mc[0]->Integral() << "\t" << cat_itr->hist_data[0]->Integral()
		<< std::endl;
    else 
      std::cout << "[DUMP NLL] " << std::setprecision(10) 
		<< cat_itr->categoryIndex1 << " " << cat_itr->categoryIndex2 
		<< "\t" << cat_itr->nll 
		<< "\t" << cat_itr->mc_events->size() << "\t" << cat_itr->data_events->size() << "\t0" 
		<< "\t" << cat_itr->hist_mc[0]->Integral() << "\t" << cat_itr->hist_data[0]->Integral()
		<< std::endl;
  }
  return;
}

