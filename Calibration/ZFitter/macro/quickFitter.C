quickFitter(TString file){//This is the idea of the iterative fit inside macro/macrofit.c (IterMinimumFit is called)
 
  Double_t asymmetricParabola(Double_t* x,Double_t* par)
  {
    Double_t value=0.;
    if (x[0]<=par[1])
      value=par[2]*(x[0]-par[1])*(x[0]-par[1]);
    else
      value=par[3]*(x[0]-par[1])*(x[0]-par[1]);
    return value+par[0];
  }
  
  TF1* f=new TF1("f","pol2",0.98,1.02); 
  g->Fit(f,"OFR+","",0.98,1.02);                                                                                                                                            
  TF1* fun=new TF1("fun",asymmetricParabola,0.98,1.02,4);                                                                                                                   
  fun->SetParameter(0,0);                                                                                                                                                   
  double minX=f->GetMinimumX();                                                                                                                                             
  fun->SetParameter(1,minX);                                                                                                                                                
  double sigma=1./sqrt(2* f->GetParameter(2));                                                                                                                              
  fun->SetParameter(2,1/(2*sigma*sigma));                                                                                                                                   
  fun->SetParameter(3,1/(2*sigma*sigma));                                                                                                                                   
  g->Fit(fun,"OFR+","",0.98,1.02);
}
