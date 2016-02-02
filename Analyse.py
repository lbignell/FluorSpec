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
import matplotlib.pyplot as plt
import time

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
        print("Initializing FluorSpecReader at {0}".format(time.time()))
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

    def ApplyCorrFileToRaw(self, rawspec, data, key, bckgnd=0, extracorr=None, MakePlots=False):
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
        if MakePlots and extracorr is None:
            #Plot the raw data.
            plt.figure()
            plt.plot(data.WL, rawspec, label='Raw Data')
            plt.plot(data.WL, CorrData, label='Corrected using {0}'.format(key))
            plt.legend()

        if extracorr is not None:
            if extracorr.RunType.name!='Synchronous':
                print('ERROR!! The extracorr run type is not synchronous!')
                return
            CorrVals = np.interp(data.WL, extracorr.WL, extracorr.Spec)
            CorrData = np.divide(CorrData,CorrVals)
            if MakePlots:
                plt.figure()
                plt.plot(data.WL, rawspec, label='Raw Data')
                plt.plot(data.WL, CorrData, label='Corrected using {0}'.format(key))
                plt.plot(data.WL, CorrData,
                         label='Corrected with Sync Scan file {0}'.format(
                         extracorr.FilePath))
                plt.legend()
        return CorrData
    
    def CalculateQY_2MM(self, corrspec_fluor, corrspec_solvent, fluor, solvent,
                        scat_start, scat_end, em_start, em_end, use_solvent_BL=False):
        '''
        Calculate QY using 2 measurement method.

        This will automatically subtract the baseline, by fitting lines between
        the start and end of the integration ranges.
        
        Arguments:
        - corrected spectrum for fluorophore and solvent as lists, obtained with
        ApplyCorrFileToRaw.
        - the PTI_Data objects for the fluorophore and solvent
        - the integration ranges for the scatter peak and emission spectrum.
        '''
        ScatStartIdx_Fluor = fluor.WL.index(np.interp(scat_start, fluor.WL, fluor.WL))
        ScatEndIdx_Fluor = fluor.WL.index(np.interp(scat_end, fluor.WL, fluor.WL))
        ScatStartIdx_Solvent = solvent.WL.index(np.interp(scat_start, solvent.WL, solvent.WL))
        ScatEndIdx_Solvent = solvent.WL.index(np.interp(scat_end, solvent.WL, solvent.WL))
        EmStartIdx_Fluor = fluor.WL.index(np.interp(em_start, fluor.WL, fluor.WL))
        EmEndIdx_Fluor = fluor.WL.index(np.interp(em_end, fluor.WL, fluor.WL))

        Scat_BL_grad_Fluor = (np.mean(corrspec_fluor[ScatEndIdx_Fluor+1:ScatEndIdx_Fluor+6]) - \
                            np.mean(corrspec_fluor[ScatStartIdx_Fluor-6:ScatStartIdx_Fluor-1]))/ \
                            (ScatEndIdx_Fluor - ScatStartIdx_Fluor)
        Scat_BL_grad_Solvent = (np.mean(corrspec_solvent[ScatEndIdx_Solvent+1:ScatEndIdx_Solvent+6]) - \
                            np.mean(corrspec_solvent[ScatStartIdx_Solvent-6:ScatStartIdx_Solvent-1]))/ \
                            (ScatEndIdx_Solvent - ScatStartIdx_Solvent)

        Scat_BL_const_Fluor = corrspec_fluor[ScatEndIdx_Fluor] - Scat_BL_grad_Fluor*fluor.WL[ScatEndIdx_Fluor]
        Scat_BL_const_Solvent = corrspec_solvent[ScatEndIdx_Solvent] - Scat_BL_grad_Solvent*solvent.WL[ScatEndIdx_Solvent]

        Scat_BL_Fluor = np.add(np.multiply(fluor.WL,Scat_BL_grad_Fluor),Scat_BL_const_Fluor)
        Scat_BL_Solvent = np.add(np.multiply(solvent.WL,Scat_BL_grad_Solvent),Scat_BL_const_Solvent)

        if use_solvent_BL:
            #do it assuming same WL range for now.
            Em_BL_Fluor = corrspec_solvent
        else:
            Em_BL_grad_Fluor = (np.mean(corrspec_fluor[EmEndIdx_Fluor+1:EmEndIdx_Fluor+6]) - \
                                np.mean(corrspec_fluor[EmStartIdx_Fluor-6:EmStartIdx_Fluor-1]))/ \
                                (EmEndIdx_Fluor - EmStartIdx_Fluor)
            Em_BL_const_Fluor = corrspec_fluor[ScatEndIdx_Fluor] - Em_BL_grad_Fluor*fluor.WL[ScatEndIdx_Fluor]
            Em_BL_Fluor = np.add(np.multiply(fluor.WL,Em_BL_grad_Fluor),Em_BL_const_Fluor)

        N_emitted = sum(np.subtract(corrspec_fluor[EmStartIdx_Fluor:EmEndIdx_Fluor],
                                    Em_BL_Fluor[EmStartIdx_Fluor:EmEndIdx_Fluor]))
        N_Tot_empty = sum(np.subtract(corrspec_solvent[ScatStartIdx_Solvent:ScatEndIdx_Solvent],
                                      Scat_BL_Solvent[ScatStartIdx_Solvent:ScatEndIdx_Solvent]))
        N_Tot_sample = sum(np.subtract(corrspec_fluor[ScatStartIdx_Fluor:ScatEndIdx_Fluor],
                                       Scat_BL_Fluor[ScatStartIdx_Fluor:ScatEndIdx_Fluor]))
        QY = N_emitted/(N_Tot_empty - N_Tot_sample)

        print("# emitted = {0}, # tot (no sample) = {1}, # tot (sample) = {2}, QY = {3}".format(
        N_emitted, N_Tot_empty, N_Tot_sample, QY))
        plt.figure()        
        plt.plot(fluor.WL, corrspec_fluor, 'b', label='fluor spec')
        plt.plot(solvent.WL, corrspec_solvent, 'r', label='solvent spec')
        plt.plot(fluor.WL, Scat_BL_Fluor, 'g', label='fluor scattering baseline')
        plt.plot(fluor.WL, Em_BL_Fluor, 'c', label='fluor emission baseline')
        plt.plot(solvent.WL, Scat_BL_Solvent, 'm', label='solvent scattering baseline')
        plt.legend()
        return
