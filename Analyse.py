#!/usr/bin/env python2
'''
Code for handling the PTI Fluorescence Spectrometer data analysis.
Some of this code borrows from David's github repository: https://github.com/davidjaffe/QY
'''

import sys
import os
import numpy as np
from ROOT import TH1D, TFile, gROOT, TCanvas, TGraph, TLegend, TGraphErrors, gStyle
from enum import Enum
import PTI_Data

class FluorSpecReader():
    '''
    Read and process data from the PTI Flourescence Spectrometer.
    
    It is assumed that the data is a text file generated either by 'Export Session',
    or 'Export Trace'.
    '''
    CorrFiles = {'emcorri':'correction_data\\emcorri.txt',
                 'emcorr-sphere':'correction_data\\emcorr-sphere.txt',
                 'emcorr-sphere-quanta':'correction_data\\emcorr-sphere-quanta.txt',
                 'excorr':'correction_data\\excorr.txt'}    
    def __init__(self):
        #Initialize parameters
        print("Initializing!")
    def GetCorrData(self, key):
        '''
        Return spectral correction data object.
        
        key may be: 'emcorri', 'emcorr-sphere', 'emcorr-sphere-quanta', or 'excorr'.
        '''
        return PTI_Data.PTI_Data(self.CorrFiles[key])
    