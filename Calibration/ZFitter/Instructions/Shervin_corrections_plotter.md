```
root -l -b tmp/d_chain.root macro/load_singleFile.C 
TString commonCut="(((((eleID[0] & 2)==2)&&((eleID[1] & 2)==2))&&((energySCEle_corr[0]/cosh(etaSCEle[0]) >= 25)&&(energySCEle_corr[1]/cosh(etaSCEle[1]) >= 25)))&&(HLTfire==1))&&(recoFlagsEle[0] > 1 && recoFlagsEle[1] > 1)";
data->Draw("scaleEle>>h","(smearerCat[0]==0)&&"+commonCut);
h->GetXaxis()->SetTitle("scale EB High mass")
c1->SaveAs("~/scratch1/www/Pt1Pt2/Likelihood/ptRatio_pt2Sum_withoutCorr/scale_EB_highMass.png")
data->Draw("scaleEle>>h_EE","(smearerCat[0]==2)&&"+commonCut);
h_EE->GetXaxis()->SetTitle("scale EE High mass")
c1->SaveAs("~/scratch1/www/Pt1Pt2/Likelihood/ptRatio_pt2Sum_withoutCorr/scale_EE_highMass.png")
h->GetMean()
h_EE->GetMean()
```
```
root -l -b tmp/s_chain.root macro/load_singleFile.C 
TString commonCut="(((((eleID[0] & 2)==2)&&((eleID[1] & 2)==2))&&((energySCEle_corr[0]/cosh(etaSCEle[0]) >= 25)&&(energySCEle_corr[1]/cosh(etaSCEle[1]) >= 25)))&&(HLTfire==1))&&(recoFlagsEle[0] > 1 && recoFlagsEle[1] > 1)";
data->Draw("smearEle>>h","(smearerCat[0]==0)&&"+commonCut);
h->GetXaxis()->SetTitle("sigma EB High mass")
c1->SaveAs("~/scratch1/www/Pt1Pt2/Likelihood/ptRatio_pt2Sum_withoutCorr/sigma_EB_highMass.png")
data->Draw("smearEle>>h_EE","(smearerCat[0]==2)&&"+commonCut);
h_EE->GetXaxis()->SetTitle("sigma EE High mass")
c1->SaveAs("~/scratch1/www/Pt1Pt2/Likelihood/ptRatio_pt2Sum_withoutCorr/sigma_EE_highMass.png")
h->GetRMS()
h_EE->GetRMS()
```