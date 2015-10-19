#ifndef roosmearer_hh
#define roosmearer_hh

#include <iomanip>

#include <TChain.h>
#include <vector>
#include <TH1F.h>
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>
#include <Math/Util.h>
#include <TStopwatch.h>

#include <RooAbsReal.h>
#include <RooSetProxy.h>
#include <RooStats/MarkovChain.h>
#include "SmearingImporter.hh"

class ZeeCategory {
  private:
  int Variables_=-1;//How many variables per cat: default is -1 to force you to settle it correctly

public:

  int GetNumberVariables(){
    return Variables_;
  }

  void SetNumberVariables(int variables){
    //Once you know the number of variables you are dealing with you can initialize your vector of histograms
    //This is called by RooSmearer::InitCategories
    Variables_=variables;
    for(int var=0;(unsigned)var<GetNumberVariables();var++){
      hist_data.push_back(NULL);
      smearHist_data.push_back(NULL);
      hist_mc.push_back(NULL);
      smearHist_mc.push_back(NULL);
    }
    if(hist_data.size()!=GetNumberVariables()){
      std::cout<<"[WARNING] ZeeCategory is not allocated correctly in its size. Have a look at interface/RooSmearer.hh, in which ZeeCategory is defined"<<std::endl;
    }
  }

  inline  ZeeCategory(){};
  
  inline ~ZeeCategory(){
    
    data_events=NULL;
    mc_events=NULL;
  };


public:
  //! Members of ZeeCategory
  zee_events_t *data_events;
  zee_events_t *mc_events;

  int categoryIndex1, categoryIndex2;
  TString categoryName1;
  TString categoryName2;
  RooArgSet pars1;
  RooArgSet pars2;

  // old values
  double scale1,constant1, alpha1;
  double scale2,constant2, alpha2;

  std::vector<int> nBins;
  std::vector<float> targetVariable_min;
  std::vector<float> targetVariable_max;
  
  //TH1F * hist_data;
  //TH1F * smearHist_data;
  //TH1F *hist_mc;
  //TH1F *smearHist_mc;

  std::vector<TH1F *> hist_data;//(Variables,NULL);
  std::vector<TH1F *> smearHist_data;//(Variables,NULL);
  std::vector<TH1F *> hist_mc;//(Variables,NULL);
  std::vector<TH1F *> smearHist_mc;//(Variables,NULL);

  bool active;
  unsigned int nSmearToy;
  unsigned int nLLtoy;

  double nll, nllRMS;
};

// parameters params are provided externally to leave flexibiility in
// their definition and dependencies (can be  RooRealVar or
// RooFormulaVar)

class RooSmearer: public RooAbsReal {

public:
  ~RooSmearer(void);

  inline RooSmearer(){};
  RooSmearer(const RooSmearer& c, const char *newname);
  virtual TObject* clone(const char* newname) const
   {
     return new RooSmearer(*this,newname);
   };

  Double_t evaluate(void) const;

  RooSmearer(const char *name, TChain *data_chain_, 
	     TChain *signal_chain_,
	     TChain *bkg_chain_, 
	     std::vector<TString> regionList,
	     std::vector<RooArgSet > params,
	     RooArgSet parset,
	     TString energyBranchName="energySCEle_regrCorr_ele"
	     );


  // Setting options
  void AutoNSmear(ZeeCategory& category);
  void AutoNBins(ZeeCategory& category);

  // interface functions to SmearingImporter
  inline void SetPuWeight(bool usePuReweight){importer.SetPuWeight(usePuReweight);};
  inline void SetR9Weight(bool useR9Reweight){importer.SetR9Weight(useR9Reweight);};
  inline void SetPtWeight(bool usePtReweight){importer.SetPtWeight(usePtReweight);};
  inline void SetFsrWeight(bool value){ importer.SetFsrWeight(value); };
  inline void SetWeakWeight(bool value){importer.SetWeakWeight(value);};
  inline void SetZPtWeight(bool useZPtReweight){importer.SetZPtWeight(useZPtReweight);};
  inline void SetOnlyDiagonal(bool value){importer.SetOnlyDiagonal(value);};
  inline void SetSmearingEt(bool value){importer.SetSmearingEt(value);};
  inline void SetPdfSystWeight(int value){importer.SetPdfSystWeight(value);};

  void SetNSmear(unsigned int n_smear=0, unsigned int nlltoy=0);

  inline void SetToyScale(float scaleToy=1.01, float constTermToy=0.01){
    importer._scaleToy=scaleToy;
    importer._constTermToy=constTermToy;
  }
  inline void SetEleID(TString value){importer.SetEleID(value);};
  inline void SetCommonCut(TString cut){importer.SetCommonCut(cut);};

  int GetNumberVariables(){return (int)targetVariable_.size();}  
  inline void SetTargetVariable(std::vector<TString> targetVariable,TString configuration){
    targetVariable_=targetVariable;//
    importer.SetTargetVariable(targetVariable);
    importer.SetConfiguration(configuration);//This sets if 1/2 are leading/subleading or random
  }

  inline void SetHistBinning(std::vector<float> min, std::vector<float> max, std::vector<float> width){
    //! Those are attributes of RooSmearer
    targetVariable_min_=min;
    targetVariable_max_=max;
    targetVariable_bin_=width;
    for(int var=0;(unsigned)var < targetVariable_min_.size();var++){
      nBins_.push_back((int) ((targetVariable_max_[var] - targetVariable_min_[var])/targetVariable_bin_[var]));//int(3)=2 (?)
    }
    return;
  }

  inline std::vector<int> GetHistBinning(){return nBins_;}
  /// Initialize the categories: import from the tree
  void Init(TString commonCut, TString eleID, Long64_t nEvents=0, bool mcToy=false, bool externToy=true,TString initFile="");

  inline RooDataSet *GetMarkovChainAsDataSet(){
    return _markov.GetAsDataSet();
  };
  inline RooDataSet *SetDataSet(TString name="profile", TString title="", double nllMin_=0){
    if(dataset!=NULL){
      std::cerr << "[WARNING] Removing last dataset: " << std::endl;
      dataset->Print();
      delete dataset;
    }
    dataset = new RooDataSet(name, title, RooArgSet(_paramSet,nllVar));
    return dataset;
  };
  inline RooDataSet *GetDataSet(void){
    return dataset;
  }


private:
  //! private members of RooSmearer
  TChain *_data_chain, *_signal_chain;
  SmearingImporter importer;
  std::vector<zee_events_t> mc_events_cache;
  std::vector<zee_events_t> data_events_cache;
  std::vector<RooArgSet> _params_vec; //, _truth_params_vec;

  RooSetProxy _paramSet;
  RooArgSet *truthSet, pullArgs;

  std::vector<float> targetVariable_min_;
  std::vector<float> targetVariable_max_;
  std::vector<float> targetVariable_bin_;
  std::vector<int> nBins_;
  std::vector<TString> targetVariable_;

public:
  float deltaNLLMaxSmearToy;
  unsigned int _deactive_minEventsDiag;
  unsigned int _deactive_minEventsOffDiag;
  double nllMin;
  unsigned int _nSmearToy;

private:

  unsigned int _nLLtoy;
  TRandom3* rgen_;
  TStopwatch *myClock;
  
  double lastNLL;
  double lastNLLrms;
  double nllBase;
  RooRealVar nllVar;

public:
  bool _isDataSmeared;
  bool _autoBin;
  bool _autoNsmear;
  bool smearscan;
  RooDataSet *dataset;

  RooStats::MarkovChain _markov;

private:
  void SetCache(Long64_t nEvents=0, bool cacheToy=false, bool externToy=true);
  void InitCategories(bool mcToy=false);

  double smearedEnergy(double *smear, unsigned int nGen, float ene,float scale,float alpha,float constant, const float *fixedSmearings=NULL) const;
  void SetSmearedHisto(const zee_events_t& cache, 
		       RooArgSet pars1, RooArgSet pars2, 
		       TString categoryName1, TString categoryName2, unsigned int nSmearToy,
		       std::vector<TH1F *> hist) const;

  void SetHisto(const zee_events_t& cache, std::vector<TH1F *> hist) const;
  //void SetAutoBin(ZeeCategory& category, double min, double max); // set using statistics
  void ResetBinning(ZeeCategory& category);
  bool isCategoryChanged(ZeeCategory& category, bool updateVar=true) const;
  

  double getLogLikelihood(std::vector<TH1F*> data, std::vector<TH1F*> prob) const;
  void UpdateCategoryNLL(ZeeCategory& cat, unsigned int nLLtoy, bool multiSmearToy=true);
  

  int Trag_eq(int row, int col, int N) const;

public:
  std::vector<ZeeCategory> ZeeCategories;
  typedef std::vector<ZeeCategory> ZeeCategoryCollection;

  std::vector<TH1F *> GetSmearedHisto(ZeeCategory& category, bool isMC,
			bool smearEnergy, bool forceNew=false, bool multiSmearToy=true);

  
  double getCompatibility(bool forceUpdate=false) const;
  void DumpNLL(void) const;
  inline RooArgSet GetParams(void){return _paramSet;};

  inline double GetNllRMS(){return lastNLLrms;};
};

#endif
