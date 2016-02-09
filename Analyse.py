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
import math

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
        print("Initializing FluorSpecReader at {0}".format(time.asctime(time.localtime())))
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

    def ApplyCorrFileToRaw(self, data, key, bckgnd=0, extracorr=None, MakePlots=False, factor=None):
        '''
        Take raw data as input and return the corrected spectrum.
        
        Arguments:
        Raw spectrum as list.
        PTI_Data object for data.
        key may be 'emcorri', 'emcorr-sphere', 'emcorr-sphere-quanta', or 'excorr'.
        bckgnd is the background (list) to be subtracted from the raw data.
        extracorr is an optional argument to apply an extra correction (to be divided)
        using synchronous scan data in the form of a PTI_Data object.
        Be sure that the right correction that was applied to the synchronous scan is
        the same as the key argument.
        *factor* optionally applies a multiplicative factor to the data (for comparing data with different 
        integration times or sensitivities).
        '''
        if data.FileType.value > 1:
            rawspec = data.Trace
            uspec = data.UTrace
        elif data.FileType.value == 1:
            rawspec = data.SpecRaw
            uspec = data.USpecRaw
        else:
            print("Analyse.ApplyCorrFileToRaw ERROR!! Bad file")
            return
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
        else:
            CorrVals = [1 for i in len(rawspec)]
            
        CorrData = np.multiply(np.subtract(rawspec, bckgnd),
                               CorrVals)
        Ubckgnd = np.sqrt(bckgnd)
        UCorrData = np.multiply(np.sqrt(np.add(np.power(uspec,2),
                                               np.power(Ubckgnd,2))),
                                CorrVals)
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
            #The mess below corrects for the measured difference between the sphere
            #response and the file correction, as well as its uncertainty.
            extracorr_vals = np.interp(data.WL, extracorr.WL, extracorr.Spec)
            CorrData_extracorr = np.divide(CorrData,extracorr_vals)
            Raw_extracorr = np.interp(data.WL, extracorr.WL, extracorr.SpecRaw)
            URaw_extracorr = np.sqrt(Raw_extracorr)
            Uextracorr = np.multiply(extracorr_vals,np.divide(URaw_extracorr, Raw_extracorr))
            UCorrData = np.sqrt(np.add(np.power(np.multiply(np.divide(1,extracorr_vals), UCorrData),2),
                               np.power(np.multiply(np.divide(CorrData_extracorr,extracorr_vals), Uextracorr),2)))
            CorrData = CorrData_extracorr
            if MakePlots:
                plt.figure()
                plt.plot(data.WL, rawspec, label='Raw Data')
                plt.plot(data.WL, CorrData, label='Corrected using {0}'.format(key))
                plt.plot(data.WL, CorrData,
                         label='Corrected with Sync Scan file {0}'.format(
                         extracorr.FilePath))
                plt.legend()
        if factor is not None:
            CorrData = np.multiply(CorrData, factor)
            UCorrData = np.multiply(UCorrData, factor)
        data.RegisterCorrSpec(CorrData,UCorrData)
        return CorrData, UCorrData
    
    def CalculateQY_2MM(self, fluor, solvent, scat_start, scat_end, em_start, em_end, use_solvent_BL=False,
                        dilute=None, normWL=None, verbose=False, avglen=4):
        '''
        Calculate QY using 2 measurement method.

        This will automatically subtract the baseline, by fitting lines between
        the start and end of the integration ranges.
        
        Arguments:
        - the PTI_Data objects for the fluorophore and solvent. Note that the correction has to be registered.
        - the integration ranges for the scatter peak and emission spectrum.
        - optionally use the solvent spectrum for the fluorescence baseline corr.
        - corrected spectrum and PTI_Data object for a dilute fluor spec (for reabsorption correction).
        - flag to print info to console and create a plot.
        '''
        ScatStartIdx_Solvent = solvent.WL.index(np.interp(scat_start, solvent.WL, solvent.WL))
        ScatEndIdx_Solvent = solvent.WL.index(np.interp(scat_end, solvent.WL, solvent.WL))
        ScatStartIdx_Fluor = fluor.WL.index(np.interp(scat_start, fluor.WL, fluor.WL))
        ScatEndIdx_Fluor = fluor.WL.index(np.interp(scat_end, fluor.WL, fluor.WL))
        EmStartIdx_Fluor = fluor.WL.index(np.interp(em_start, fluor.WL, fluor.WL))
        EmEndIdx_Fluor = fluor.WL.index(np.interp(em_end, fluor.WL, fluor.WL))
        #Calculate baselines
        Scat_BL_Fluor = self.CalcStraightLine(fluor.WL,
                                              fluor.SpecCorrected,
                                              ScatStartIdx_Fluor,
                                              ScatEndIdx_Fluor,
                                              avglen=avglen)
        Scat_BL_Solvent = self.CalcStraightLine(solvent.WL,
                                              solvent.SpecCorrected,
                                              ScatStartIdx_Solvent,
                                              ScatEndIdx_Solvent,
                                              avglen=avglen)
        if use_solvent_BL:
            #do it assuming same WL range for now.
            Em_BL_Fluor = solvent.SpecCorrected
        else:
            Em_BL_Fluor = self.CalcStraightLine(fluor.WL,
                                                fluor.SpecCorrected,
                                                EmStartIdx_Fluor,
                                                EmEndIdx_Fluor,
                                                avglen=avglen)

        N_emitted = sum(np.subtract(fluor.SpecCorrected[EmStartIdx_Fluor:EmEndIdx_Fluor],
                                    Em_BL_Fluor[EmStartIdx_Fluor:EmEndIdx_Fluor]))
        UN_emitted = np.sqrt(sum(np.power(fluor.USpecCorrected[EmStartIdx_Fluor:EmEndIdx_Fluor],2)))
        N_Tot_empty = sum(np.subtract(solvent.SpecCorrected[ScatStartIdx_Solvent:ScatEndIdx_Solvent],
                                      Scat_BL_Solvent[ScatStartIdx_Solvent:ScatEndIdx_Solvent]))
        UN_Tot_empty = np.sqrt(sum(np.power(solvent.USpecCorrected[ScatStartIdx_Solvent:ScatEndIdx_Solvent],2)))
        N_Tot_sample = sum(np.subtract(fluor.SpecCorrected[ScatStartIdx_Fluor:ScatEndIdx_Fluor],
                                       Scat_BL_Fluor[ScatStartIdx_Fluor:ScatEndIdx_Fluor]))
        UN_Tot_sample = np.sqrt(sum(np.power(fluor.USpecCorrected[ScatStartIdx_Fluor:ScatEndIdx_Fluor],2)))
        if((dilute is not None)) and (normWL is not None):
            if dilute.SpecCorrected is not None:
                w, Uw = self.CalcReabsProb(fluor, em_start, em_end,
                                           dilute.WL.index(np.interp(normWL, dilute.WL, dilute.WL)),
                                           Em_BL_Fluor, dilute, verbose)
        else:
            w = 0
            Uw = 0

        QY = N_emitted/(N_Tot_empty - N_Tot_sample)
        UQY = np.sqrt((1/(N_Tot_empty-N_Tot_sample)**2)*(UN_emitted**2) + 
                      ((QY/(N_Tot_empty-N_Tot_sample))**2)*(UN_Tot_empty**2 + UN_Tot_sample**2))
        QY_w = QY/(1- w + w*QY)
        UQY = np.sqrt(((QY_w/QY)**2)*(UQY**2) +
                      (((1-QY)*QY_w/(1-w+w*QY))**2)*(Uw**2))
        QY = QY_w
        if verbose:
            print("Quantum Yield: \n # emitted = {0} +/- {1}, \
            # tot (no sample) = {2} +/- {3}, \
            # tot (sample) = {4} +/- {5}, \
            QY = {6} +/- {7}".format(
                N_emitted, UN_emitted, N_Tot_empty, UN_Tot_empty, N_Tot_sample, UN_Tot_sample, QY, UQY))
            plt.figure()        
            plt.errorbar(fluor.WL, fluor.SpecCorrected, yerr=fluor.USpecCorrected, color='b', label='fluor spec')
            plt.errorbar(solvent.WL, solvent.SpecCorrected, yerr=solvent.USpecCorrected, color='r', label='solvent spec')
            plt.plot(fluor.WL, Scat_BL_Fluor, 'g', label='fluor scattering baseline')
            plt.plot(fluor.WL, Em_BL_Fluor, 'c', label='fluor emission baseline')
            plt.plot(solvent.WL, Scat_BL_Solvent, 'm', label='solvent scattering baseline')
            plt.legend()
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Fluorescence Intensity (AU)')
            plt.title('file: ' + str(fluor.FilePath.split('\\')[-1]) +
                '\n Excitation' + str(fluor.ExRange) +
                ' nm, Emission ' + str(fluor.EmRange) + ' nm')
            plt.figure()
            plt.errorbar(solvent.WL[ScatStartIdx_Solvent:ScatEndIdx_Solvent], np.subtract(solvent.SpecCorrected[ScatStartIdx_Solvent:ScatEndIdx_Solvent],
                                                                                          Scat_BL_Solvent[ScatStartIdx_Solvent:ScatEndIdx_Solvent]),
                         yerr=solvent.USpecCorrected[ScatStartIdx_Solvent:ScatEndIdx_Solvent],
                         color='r', label='solvent spec')
            plt.errorbar(fluor.WL[ScatStartIdx_Fluor:ScatEndIdx_Fluor], np.subtract(fluor.SpecCorrected[ScatStartIdx_Fluor:ScatEndIdx_Fluor],
                                                                                    Scat_BL_Fluor[ScatStartIdx_Fluor:ScatEndIdx_Fluor]),
                         yerr=fluor.USpecCorrected[ScatStartIdx_Fluor:ScatEndIdx_Fluor],
                         color='b', label='fluor spec, scatter')
            plt.errorbar(fluor.WL[EmStartIdx_Fluor:EmEndIdx_Fluor], np.subtract(fluor.SpecCorrected[EmStartIdx_Fluor:EmEndIdx_Fluor],
                                                                                Em_BL_Fluor[EmStartIdx_Fluor:EmEndIdx_Fluor]),
                         yerr=fluor.USpecCorrected[EmStartIdx_Fluor:EmEndIdx_Fluor],
                         color='g', label='fluor spec, emission')
            plt.plot((fluor.WL[ScatStartIdx_Fluor],fluor.WL[EmEndIdx_Fluor]), (0,0), 'k')
            plt.legend()
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Fluorescence Intensity (AU)')
            plt.title('file: ' + str(fluor.FilePath.split('\\')[-1]) +
                '\n Excitation' + str(fluor.ExRange) +
                ' nm, Emission ' + str(fluor.EmRange) + ' nm')
        return QY, UQY

    def CalcStraightLine(self, WL, spec, startidx, endidx, avglen=4):
        gradient = (np.mean(spec[(endidx):(endidx+avglen)]) - np.mean(spec[(startidx-avglen):(startidx)]))/(WL[endidx+int(avglen/2)] - WL[startidx-int(avglen/2)])
        #print('spec[(endidx):(endidx+6)] = {0}, spec[(startidx):(startidx-6)] = {1}, WL[endidx+3] = {2}, WL[startidx-3] = {3}'.format(
         #   spec[(endidx):(endidx+6)], spec[(startidx):(startidx-6)], WL[endidx+3], WL[startidx-3]))
        #gradient = (spec[endidx]-spec[startidx])/(WL[endidx] - WL[startidx])
        const = np.mean(spec[endidx:(endidx+avglen)]) - gradient*WL[endidx]
        #print('gradient = {0}, constant = {1}'.format(gradient, const))
        return np.add(np.multiply(WL,gradient),const)

    def CalcReabsProb(self, sphere, em_start, em_end, normWL,
                      Em_BL_Fluor, dilute, verbose=False):
        '''
        Calculate the reabsorption probability (necessary correction).
        
        Arguments:
        - Corrected fluorescence spec in sphere.
        - PTI_Data object for above.
        - Emission integration range start, then end.
        - A wavelength to normalise both spectra to (needs to be somwhere that
          still has counts but minimal reabsorption).
        - Corrected fluorescence spec for a dilute (not reabsorbed) sample.
        - PTI_Data object for above.
        - Option to print results to console.
        '''
        StartIdx_Sphere = sphere.WL.index(np.interp(em_start, sphere.WL, sphere.WL))
        EndIdx_Sphere = sphere.WL.index(np.interp(em_end, sphere.WL, sphere.WL))
        StartIdx_Dilute = dilute.WL.index(np.interp(em_start, dilute.WL, dilute.WL))
        EndIdx_Dilute = dilute.WL.index(np.interp(em_end, dilute.WL, dilute.WL))
        spherespec = np.subtract(sphere.SpecCorrected, Em_BL_Fluor)
        integ_Sphere = sum(np.divide(spherespec[StartIdx_Sphere:EndIdx_Sphere],
                                     spherespec[normWL]))
        Uinteg_sphere = np.sqrt(sum(np.power(np.divide(spherespec[StartIdx_Sphere:EndIdx_Sphere],
                                                       spherespec[normWL]),2)))
        integ_Dilute = sum(np.divide(dilute.SpecCorrected[StartIdx_Dilute:EndIdx_Dilute],
                                     dilute.SpecCorrected[normWL]))
        Uinteg_Dilute = np.sqrt(sum(np.power(np.divide(dilute.SpecCorrected[StartIdx_Sphere:EndIdx_Sphere],
                                                       dilute.SpecCorrected[normWL]),2)))
        w = 1 - (integ_Sphere/integ_Dilute)
        Uw = np.sqrt((1/(integ_Dilute**2))*(Uinteg_sphere**2) +
                    ((integ_Sphere/(integ_Dilute**2))**2)*(Uinteg_Dilute**2))
        if verbose:
            print("Reabsorption calculation:\n Sphere integral = {0}, Dilute integral = {1}, 1-w = {2}, w = {3} +/- {4}".format(
                    integ_Sphere, integ_Dilute, 1-w, w, Uw))
            plt.figure()
            plt.errorbar(sphere.WL[StartIdx_Sphere:EndIdx_Sphere], np.divide(
                spherespec[StartIdx_Sphere:EndIdx_Sphere],spherespec[normWL]),
                yerr=np.divide(sphere.USpecCorrected[StartIdx_Sphere:EndIdx_Sphere],spherespec[normWL]),
                color='g', label='QY spectrum (with reabsorption)')
            plt.errorbar(dilute.WL[StartIdx_Dilute:EndIdx_Dilute], np.divide(
                dilute.SpecCorrected[StartIdx_Dilute:EndIdx_Dilute], dilute.SpecCorrected[normWL]),
                yerr=np.divide(dilute.USpecCorrected[StartIdx_Dilute:EndIdx_Dilute],dilute.SpecCorrected[normWL]),
                color='r', label='Dilute spectrum (no reabsorption)')
            plt.legend(fontsize=12)
            plt.title('file: ' + str(sphere.FilePath.split('\\')[-1]) +
                '\n Excitation' + str(sphere.ExRange) +
                ' nm, Emission ' + str(sphere.EmRange) + ' nm')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Fluorescence Intensity (AU)')
        return w, Uw