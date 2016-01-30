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
                 'excorr':'correction_data\\excorr.txt',
                 'default':None}    

    def __init__(self):
        pass

    def GetCorrData(self, key):
        '''
        Return spectral correction data object.
        
        key may be: 'emcorri', 'emcorr-sphere', 'emcorr-sphere-quanta', or 'excorr'.
        '''
        if key not in self.CorrFiles:
            print('ERROR!! Incorrect choice of correction file.')
            return
        return PTI_Data.PTI_Data(self.CorrFiles[key])

    def ApplyCorrFileToRaw(self, rawspec, data, key, bckgnd=0, extracorr=None):
        '''
        Take raw data as input and return the corrected spectrum.
        
        Arguments:
        Raw spectrum as list.
        PTI_Data object for data.
        key may be 'emcorri', 'emcorr-sphere', 'emcorr-sphere-quanta', or 'excorr'.
        bckgnd is the background to be subtracted from the raw data.
        extracorr is an optional argument to apply an extra correction (to be divided)
        using synchronous scan data in the form of a PTI_Data object.
        Be sure that the right correction that was applied to the synchronous scan is
        the same as the key argument.
        '''
        CorrData = None
        if key is not 'default':
            corr = self.GetCorrData(key)
            if corr is None:
                print('Not correcting data.')
                return
            if data.RunType.value!=corr.RunType.value:
                print('ERROR!! The correction type doesn\'t match the data type.')
                return
            CorrVals = np.interp(data.WL, corr.WL, corr.Trace, left=0, right=0)
            CorrData = np.multiply(np.subtract(rawspec, bckgnd), CorrVals)
        else:
            CorrData = np.subtract(rawspec,bckgnd)
        if extracorr is not None:
            if extracorr.RunType.name!='Synchronous':
                print('ERROR!! The extracorr run type is not synchronous!')
                return
            CorrVals = np.interp(data.WL, extracorr.WL, extracorr.Spec)
            CorrData = np.divide(CorrData,CorrVals)
        return CorrData
