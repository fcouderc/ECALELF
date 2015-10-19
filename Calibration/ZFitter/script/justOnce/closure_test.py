#! /usr/bin/python
import math
import ROOT
import numpy as np
import os
ROOT.gSystem.Load("~/rootlogon_C.so")
ROOT.rootlogon()

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
#ROOT.gStyle.SetFillStyle(ROOT.kWhite)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetFrameBorderMode(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetCanvasBorderMode(ROOT.kWhite)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadBorderMode(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetStatColor(ROOT.kWhite)
ROOT.gStyle.SetErrorX(0)

#regions=['EB']
regions=['EB','EE']

true_scales                 ={}
true_sigmas                 ={}
fitted_scales               ={}
fitted_sigmas               ={}
fitted_scales_err_down      ={}
fitted_scales_err_up        ={}
fitted_sigmas_err_down      ={}
fitted_sigmas_err_up        ={}
true_scales_check           ={}
true_sigmas_check           ={}
fitted_scales_check         ={}
fitted_sigmas_check         ={}
fitted_scales_err_down_check={}
fitted_scales_err_up_check  ={}
fitted_sigmas_err_down_check={}
fitted_sigmas_err_up_check  ={}

for region in regions:
    true_scales                 [region]=[]   
    true_sigmas                 [region]=[]   
    fitted_scales               [region]=[]   
    fitted_sigmas               [region]=[]   
    fitted_scales_err_down      [region]=[]   
    fitted_scales_err_up        [region]=[]   
    fitted_sigmas_err_down      [region]=[]   
    fitted_sigmas_err_up        [region]=[]   
    true_scales_check           [region]=[]
    true_sigmas_check           [region]=[]
    fitted_scales_check         [region]=[]
    fitted_sigmas_check         [region]=[]
    fitted_scales_err_down_check[region]=[]
    fitted_scales_err_up_check  [region]=[]
    fitted_sigmas_err_down_check[region]=[]
    fitted_sigmas_err_up_check  [region]=[]

true_scales_array                 ={}   
true_sigmas_array                 ={}   
fitted_scales_array               ={}   
fitted_sigmas_array               ={}   
fitted_scales_err_down_array      ={}   
fitted_scales_err_up_array        ={}   
fitted_sigmas_err_down_array      ={}   
fitted_sigmas_err_up_array        ={}   
true_scales_check_array           ={}
true_sigmas_check_array           ={}
fitted_scales_check_array         ={}
fitted_sigmas_check_array         ={}
fitted_scales_err_down_check_array={}
fitted_scales_err_up_check_array  ={}
fitted_sigmas_err_down_check_array={}
fitted_sigmas_err_up_check_array  ={}

scale_graph       ={}
sigma_graph       ={}
scale_check_graph ={}
sigma_check_graph ={}
canvas_scale      ={}
canvas_sigma      ={}
canvas_scale_check={}
canvas_sigma_check={}

path='/afs/cern.ch/user/g/gfasanel/new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter/test/dato/fitres/toys/scaleStep0/'

injected_sigmas=['0.00','0.005','0.01','0.0105','0.02'] #0.015!!
injected_scales=['0.98','0.985','0.99','0.995','1.00','1.005','1.01','1.015','1.02']

for region in regions:
    for injected_sigma in injected_sigmas:
        for injected_scale in injected_scales:
            with open(str(path+injected_sigma+'_0.00_'+injected_scale+'/'+region+"/outProfile_ptRatio_pt2Sum_random_"+injected_scale+'_'+injected_sigma+'_scaleStep0_Et_25_trigger_noPF-FitResult.dat')) as file_res:
                for line in file_res:  #Line is a string #split the string on whitespace, return a list
                    numbers_str = line.split() #still string
                    numbers_float=[]
                    numbers_float.append(float(numbers_str[2])) #->0: fit value
                    numbers_float.append(float(numbers_str[3])) #->1: error down
                    numbers_float.append(float(numbers_str[4])) #->2: error up
                    if not region in numbers_str[1]:
                        if "constTerm" in numbers_str[0]:                    
                            fitted_sigmas_check[region].append(numbers_float[0])
                            fitted_sigmas_err_down_check[region].append(numbers_float[1])
                            fitted_sigmas_err_up_check[region].append(numbers_float[2])
                            true_sigmas_check[region].append(0.0)
                        elif "scale" in numbers_str[0]:
                            fitted_scales_check[region].append(numbers_float[0])
                            fitted_scales_err_down_check[region].append(numbers_float[1])
                            #fitted_scales_err_up_check[region].append(numbers_float[2])
                            fitted_scales_err_up_check[region].append(0.)
                            true_scales_check[region].append(1.)
                    else:
                        if "constTerm" in numbers_str[0]:                    
                            fitted_sigmas[region].append(numbers_float[0])
                            fitted_sigmas_err_down[region].append(numbers_float[1])
                            fitted_sigmas_err_up[region].append(numbers_float[2])
                            true_sigmas[region].append(float(injected_sigma))
                        elif "scale" in numbers_str[0]:
                            fitted_scales[region].append(numbers_float[0])
                            fitted_scales_err_down[region].append(numbers_float[1])
                            fitted_scales_err_up[region].append(numbers_float[2])
                            true_scales[region].append(float(injected_scale))
    #At this point lists are full of points
    #usage of array for TGraph, otherwise it doesn't work
    true_scales_array                 [region]=np.asarray(true_scales                 [region])   
    true_sigmas_array                 [region]=np.asarray(true_sigmas                 [region])   
    fitted_scales_array               [region]=np.asarray(fitted_scales               [region])   
    fitted_sigmas_array               [region]=np.asarray(fitted_sigmas               [region])   
    fitted_scales_err_down_array      [region]=np.asarray(fitted_scales_err_down      [region])   
    fitted_scales_err_up_array        [region]=np.asarray(fitted_scales_err_up        [region])   
    fitted_sigmas_err_down_array      [region]=np.asarray(fitted_sigmas_err_down      [region])   
    fitted_sigmas_err_up_array        [region]=np.asarray(fitted_sigmas_err_up        [region])
    true_scales_check_array           [region]=np.asarray(true_scales_check           [region])
    true_sigmas_check_array           [region]=np.asarray(true_sigmas_check           [region])
    fitted_scales_check_array         [region]=np.asarray(fitted_scales_check         [region])
    fitted_sigmas_check_array         [region]=np.asarray(fitted_sigmas_check         [region])
    fitted_scales_err_down_check_array[region]=np.asarray(fitted_scales_err_down_check[region])
    fitted_scales_err_up_check_array  [region]=np.asarray(fitted_scales_err_up_check  [region])
    fitted_sigmas_err_down_check_array[region]=np.asarray(fitted_sigmas_err_down_check[region])
    fitted_sigmas_err_up_check_array  [region]=np.asarray(fitted_sigmas_err_up_check  [region])

    #scale_graph[region]=ROOT.TGraph(len(true_scales_array[region]),true_scales_array[region],fitted_scales_array[region])
    #sigma_graph[region]=ROOT.TGraph(len(true_sigmas_array[region]),true_sigmas_array[region],fitted_sigmas_array[region])
    #scale_check_graph[region]=ROOT.TGraph(len(true_scales_check_array[region]),true_scales_check_array[region],fitted_scales_check_array[region])
    #sigma_check_graph[region]=ROOT.TGraph(len(true_sigmas_check_array[region]),true_sigmas_check_array[region],fitted_sigmas_check_array[region])
    #True scales and sigmas have no errors
    dummy_scale=[]
    dummy_sigma=[]
    for i in range(0,len(true_scales_array[region])):
        dummy_scale.append(0.)
    for i in range(0,len(true_sigmas_array[region])):
        dummy_sigma.append(0.)

    scale_graph[region]=ROOT.TGraphAsymmErrors(len(true_scales_array[region]),true_scales_array[region],fitted_scales_array[region],np.asarray(dummy_scale),np.asarray(dummy_sigma),fitted_scales_err_down_array[region],fitted_scales_err_up_array[region])
    sigma_graph[region]=ROOT.TGraphAsymmErrors(len(true_sigmas_array[region]),true_sigmas_array[region],fitted_sigmas_array[region],np.asarray(dummy_sigma),np.asarray(dummy_sigma),fitted_sigmas_err_down_array[region],fitted_sigmas_err_up_array[region])
    scale_check_graph[region]=ROOT.TGraphAsymmErrors(len(true_scales_check_array[region]),true_scales_check_array[region],fitted_scales_check_array[region],np.asarray(dummy_scale),np.asarray(dummy_scale),fitted_scales_err_down_check_array[region],fitted_scales_err_up_check_array[region])
    sigma_check_graph[region]=ROOT.TGraphAsymmErrors(len(true_sigmas_check_array[region]),true_sigmas_check_array[region],fitted_sigmas_check_array[region],np.asarray(dummy_sigma),np.asarray(dummy_sigma),fitted_sigmas_err_down_check_array[region],fitted_sigmas_err_up_check_array[region])

    #Plot everything
    canvas_scale[region]=ROOT.TCanvas(str("fitted_scale_"+region),str("fitted_scale_"+region))
    scale_graph[region].Draw("AP")
    scale_graph[region].GetXaxis().SetRangeUser(0.97, 1.03)
    scale_graph[region].GetYaxis().SetRangeUser(0.97, 1.03)
    scale_graph[region].SetMarkerStyle(21)
    scale_graph[region].SetMarkerColor(4)
    scale_graph[region].SetLineColor(4)
    scale_graph[region].GetXaxis().SetTitle("true scale")
    scale_graph[region].GetYaxis().SetTitle("fitted scale")
    scale_graph[region].Draw("AP")
    canvas_scale[region].Update();
    canvas_scale[region].SaveAs(str(path+region+"/fitted_scale_"+region+".png"))

    canvas_sigma[region]=ROOT.TCanvas(str("fitted_sigma_"+region),str("fitted_sigma_"+region))
    sigma_graph[region].Draw("AP")
    sigma_graph[region].GetXaxis().SetRangeUser(0.00, 0.03)
    sigma_graph[region].GetYaxis().SetRangeUser(0.00, 0.03)
    sigma_graph[region].SetMarkerStyle(21)
    sigma_graph[region].SetMarkerColor(ROOT.kRed)
    sigma_graph[region].SetLineColor(ROOT.kRed)
    sigma_graph[region].GetXaxis().SetTitle("true sigma")
    sigma_graph[region].GetYaxis().SetTitle("fitted sigma")
    sigma_graph[region].Draw("AP")
    canvas_sigma[region].Update();
    canvas_sigma[region].SaveAs(str(path+region+"/fitted_sigma_"+region+".png"))
    #check
    canvas_scale_check[region]=ROOT.TCanvas(str("fitted_scale_"+region),str("fitted_scale_"+region))
    scale_check_graph[region].Draw("AP")
    scale_check_graph[region].GetXaxis().SetRangeUser(0.97, 1.03)
    scale_check_graph[region].GetYaxis().SetRangeUser(0.97, 1.03)
    scale_check_graph[region].SetMarkerStyle(21)
    scale_check_graph[region].SetMarkerColor(4)
    scale_check_graph[region].SetLineColor(4)
    scale_check_graph[region].GetXaxis().SetTitle("true scale")
    scale_check_graph[region].GetYaxis().SetTitle("fitted scale")
    scale_check_graph[region].Draw("AP")
    canvas_scale_check[region].Update();
    canvas_scale_check[region].SaveAs(str(path+region+"/check_scale_"+region+".png"))
    canvas_sigma_check[region]=ROOT.TCanvas(str("fitted_sigma_"+region),str("fitted_sigma_"+region))
    sigma_check_graph[region].Draw("AP")
    sigma_check_graph[region].GetXaxis().SetRangeUser(0.00, 0.002)
    sigma_check_graph[region].GetYaxis().SetRangeUser(0.00, 0.03)
    sigma_check_graph[region].SetMarkerStyle(21)
    sigma_check_graph[region].SetMarkerColor(ROOT.kRed)
    sigma_check_graph[region].SetLineColor(ROOT.kRed)
    sigma_check_graph[region].GetXaxis().SetTitle("true sigma")
    sigma_check_graph[region].GetYaxis().SetTitle("fitted sigma")
    sigma_check_graph[region].Draw("AP")
    canvas_sigma_check[region].Update();
    canvas_sigma_check[region].SaveAs(str(path+region+"/check_sigma_"+region+".png"))

os.system("cp -r /afs/cern.ch/user/g/gfasanel/new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter/test/dato/fitres/toys/scaleStep0/EB/*.png ~/scratch1/www/Pt1Pt2/closure_test/EB/")
os.system("cp -r /afs/cern.ch/user/g/gfasanel/new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter/test/dato/fitres/toys/scaleStep0/EE/*.png ~/scratch1/www/Pt1Pt2/closure_test/EE/")
