#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import subprocess
from glob import glob
from collections import defaultdict
from collections import OrderedDict
from array import array

from ROOT import *

def setParameters(function,Norm,Peak):

    ## Normalisation
    function.SetParameter(0,Norm)
    function.SetParLimits(0,0.,Norm*1000.)
    
    ## 1274 KeV compton 
    function.SetParameter(1,0.4)
    function.SetParLimits(1,0.1,1)
    function.SetParameter(2,1.8*Peak)
    function.SetParLimits(2,1.7*Peak,1.9*Peak)

    ## 1274 KeV "photo-electric" 
    function.SetParameter(3,1.)
    function.SetParLimits(3,0.2,5.)
    function.SetParameter(4,2.1*Peak)
    function.SetParLimits(4,2.0*Peak,2.2*Peak)
    function.SetParameter(5,0.1*Peak)
    function.SetParLimits(5,0.03*Peak,0.2*Peak)
    function.SetParameter(17,0.4)
    function.SetParLimits(17,0.1,2.)
    function.SetParameter(18,2.)
    function.SetParLimits(18,0.1,10.)
    
    ## 511 KeV compton
    function.SetParameter(6,20.)
    function.SetParLimits(6,0.,1000.)
    function.SetParameter(7,0.4)
    function.SetParLimits(7,0.1,1)
    function.SetParameter(8,0.6*Peak)
    function.SetParLimits(8,0.05*Peak,1.3*Peak)
    
    ## Trigger turn on (Compton+BS)
    function.SetParameter(12,5.)
    function.SetParLimits(12,0,10.)
    function.SetParameter(13,5.)
    function.SetParLimits(13,0.1,10.)
    
    ## 511 KeV photoelectric
    function.SetParameter(9,15.)
    function.SetParLimits(9,0.,1000.)
    function.SetParameter(10,Peak)
    function.SetParLimits(10,0.9*Peak,1.1*Peak)
    function.SetParameter(11,0.05*Peak)
    function.SetParLimits(11,0.02*Peak,0.2*Peak)
    
    ## Back scatter peak
    function.SetParameter(14,0.01)
    function.SetParLimits(14,0.,10.)
    function.SetParameter(15,0.45*Peak)
    function.SetParLimits(15,0.3*Peak,0.6*Peak)
    function.SetParameter(16,0.07*Peak)
    function.SetParLimits(16,0.04*Peak,0.13*Peak)
    ##

def totalFunction(x,par):
    
    t = (x[0]-par[4])/par[5]
    if par[17]<0: 
        t = -t

    absAlpha = abs(par[17])
    if( t >= - absAlpha ):
        crystalball = par[3]*TMath.Exp(-0.5*t*t)
    else:
        a = TMath.Power(par[18]/absAlpha,par[18])*TMath.Exp(-0.5*absAlpha*absAlpha)
        b = par[18]/absAlpha - absAlpha
        crystalball = par[3]*(a/TMath.Power(b-t,par[18]))

    return par[0]*(1./(1+TMath.Exp(par[1]*(x[0]-par[2])))+crystalball+( 1 + TMath.Erf((x[0]-par[12])/(par[13]*TMath.Sqrt(x[0]))) )*( par[6]/(1+TMath.Exp(par[7]*(x[0]-par[8]))) + par[14]*TMath.Gaus(x[0],par[15],par[16]) )+par[9]*TMath.Gaus(x[0],par[10],par[11]))


def f_1274keV_compton(x,par):
    return par[0]*(1./(1+TMath.Exp(par[1]*(x[0]-par[2]))))

def f_1274keV_peak(x,par):
    t = (x[0]-par[1])/par[2]
    if par[3]<0: 
        t = -t

    absAlpha = abs(par[3])
    if( t >= - absAlpha ):
        crystalball = par[0]*TMath.Exp(-0.5*t*t)
    else:
        a = TMath.Power(par[4]/absAlpha,par[4])*TMath.Exp(-0.5*absAlpha*absAlpha)
        b = par[4]/absAlpha - absAlpha
        crystalball = par[0]*(a/TMath.Power(b-t,par[4]))

    return crystalball

def f_511keV_compton_times_turnon(x,par):
    return ( 1 + TMath.Erf((x[0]-par[3])/(par[4]*TMath.Sqrt(x[0]))) ) * par[0]/(1+TMath.Exp(par[1]*(x[0]-par[2])))

def f_backscatter_peak_times_turnon(x,par):
    return ( 1 + TMath.Erf((x[0]-par[3])/(par[4]*TMath.Sqrt(x[0]))) ) * par[0]*TMath.Gaus(x[0],par[1],par[2])

def f_511keV_peak(x,par):
    return par[0]*TMath.Gaus(x[0],par[1],par[2])

def f_background(x,par):
    return par[0]*(1./(1+TMath.Exp(par[1]*(x[0]-par[2])))) + ( 1 + TMath.Erf((x[0]-par[6])/(par[7]*TMath.Sqrt(x[0]))) ) * par[3]/(1+TMath.Exp(par[4]*(x[0]-par[5])))
    
def fitSpectrum(histo,function,xmin,xmax,canvas,fitres,label):

    histo.GetXaxis().SetRange(25,1000)
    peak=histo.GetBinCenter(histo.GetMaximumBin())
    norm=float(histo.GetEntries())/float(histo.GetNbinsX())
    histo.GetXaxis().SetRangeUser(xmin,xmax)
    setParameters(function,norm,peak)
    
    canvas.cd()
    histo.Draw("PE")
    goodChi2 = 0.
    previousChi2overNdf = -99.
    while goodChi2==0.:
        histo.Fit(function.GetName(),"LR+0N","",xmin,min(peak*2.4,xmax))
        print function.GetChisquare(), function.GetNDF(), function.GetChisquare()/function.GetNDF()
        if abs(function.GetChisquare()/function.GetNDF()-previousChi2overNdf)<0.01*previousChi2overNdf:
            histo.Fit(function.GetName(),"LR+","",xmin,min(peak*2.4,xmax))
            canvas.Update()
            goodChi2 = 1.
        previousChi2overNdf = function.GetChisquare()/function.GetNDF()
    print function.GetChisquare(), function.GetNDF(), function.GetChisquare()/function.GetNDF()

    fitres[(label,"peak1","mean","value")]=function.GetParameter(10)
    fitres[(label,"peak1","mean","sigma")]=function.GetParError(10)
    fitres[(label,"peak1","sigma","value")]=function.GetParameter(11)
    fitres[(label,"peak1","sigma","sigma")]=function.GetParError(11)
    fitres[(label,"peak2","mean","value")]=function.GetParameter(4)
    fitres[(label,"peak2","mean","sigma")]=function.GetParError(4)
    fitres[(label,"peak2","sigma","value")]=function.GetParameter(5)
    fitres[(label,"peak2","sigma","sigma")]=function.GetParError(5)
    fitres[(label,"backpeak","mean","value")]=function.GetParameter(15)
    fitres[(label,"backpeak","mean","sigma")]=function.GetParError(15)

    canvas.Write()

def fitSaturation(function,xmin,xmax,canvas,fitres,label):    

    canvas.cd()
    n_points = 2
    a_true_energy = [511.,1275.]
    a_err_true_energy = [0.,0.]
    a_meas_energy = [fitResults[(label,"peak1","mean","value")],
                     fitResults[(label,"peak2","mean","value")]]
    a_err_meas_energy = [fitResults[(label,"peak1","mean","sigma")],
                         fitResults[(label,"peak2","mean","sigma")]]
    g_sat = TGraphErrors(n_points,
                         array('d',a_true_energy),array('d',a_meas_energy),
                         array('d',a_err_true_energy),array('d',a_err_meas_energy))
    g_sat.GetXaxis().SetLimits(xmin,xmax)
    g_sat.GetYaxis().SetLimits(0.,fitResults[(label,"peak2","mean","value")]*1.2)
    g_sat.SetMarkerStyle(20)   
 
    function.SetParameter(0,fitResults[(label,"peak2","mean","value")])
    function.SetParLimits(0,fitResults[(label,"peak2","mean","value")]*0.9,fitResults[(label,"peak2","mean","value")]*3)
    function.SetParameter(1,0.0007)
    function.SetParLimits(1,0.00001,0.01)
    
    g_sat.Draw("APE")
    g_sat.Fit(function.GetName(),"WR","",xmin,xmax)
    g_sat.GetXaxis().UnZoom()
    g_sat.GetYaxis().UnZoom()

    fitres[(label,"peak12","alpha","value")]=function.GetParameter(0)
    fitres[(label,"peak12","alpha","sigma")]=function.GetParError(0)
    fitres[(label,"peak12","beta","value")]=function.GetParameter(1)
    fitres[(label,"peak12","beta","sigma")]=function.GetParError(1)

    canvas.Write()

###############################################
gROOT.SetBatch(True)

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetOptFit(1111111)
gStyle.SetStatH(0.09)


################################################
## 6) Fit energy spectra
################################################

## Setup
minEnergy = 4
maxEnergy = 90
n_paramameters = 19

fitResults = {}

f=TFile.Open("Spectrum.root","read")
h1_energy_pixel=f.Get("HistoCh59N0")
h1_energyTot_bar=f.Get("HistoCh288N0")


c1_energy_pixel = TCanvas("c1_energy_pixel", "", 900, 700)
c1_energyTot_bar = TCanvas("c1_sat_pixel", "", 900, 700)

## Pixel (ref)
fTot_pixel = TF1("fTot_pixel",totalFunction,minEnergy,maxEnergy,n_paramameters)
fTot_pixel.SetNpx(1000)
fitSpectrum(h1_energy_pixel,fTot_pixel,minEnergy,maxEnergy,c1_energy_pixel,fitResults,"pixel")

## Bar (ref)
fTot_bar = TF1("fTot_bar",totalFunction,minEnergy,maxEnergy,n_paramameters)
fTot_bar.SetNpx(1000)
fitSpectrum(h1_energyTot_bar,fTot_bar,minEnergy,maxEnergy,c1_energyTot_bar,fitResults,"bar")

outputFile = TFile("SpectrumFile.root","RECREATE");
outputFile.cd()
fTot_pixel.Write()
fTot_bar.Write()

outputFile.Save()
outputFile.Close()
