#!/usr/bin/env python2
'''
These classes handle data from the PTI spectrometer.
TEXT data are assumed (not the .gx* nonsense).
'''
import os
from enum import Enum
import time

class PTI_Data:
    '''PTI spectrometer data class.'''
    RunTypes = Enum('RunType', 'Unknown Emission Excitation Synchronous')
    FileTypes = Enum('FileType', 'Unknown Session Trace Group')
    def __init__(self, fname):
        #Get the file as an object.
        self.FilePath = fname
        with file(self.FilePath, 'r') as thefile:
            firstline = thefile.readline()
            if '<Session>' in firstline:
                self.FileType = self.FileTypes.Session
            elif '<Trace>' in firstline:
                self.FileType = self.FileTypes.Trace
            elif '<Group>' in firstline:
                self.FileType = self.FileTypes.Group
            else:
                print("ERROR!! Unknown file format.")
                self.FileType = self.FileTypes.Unknown
                self.SuccessfullyRead = False
                return

        self.SuccessfullyRead = self.ReadHeaderInfo()
        self.WL = [0]*self.NumSamples
        if self.FileType == self.FileTypes.Session:
            self.Spec = [0]*self.NumSamples
            self.SpecRaw = [0]*self.NumSamples
            self.ExCorr = [0]*self.NumSamples #Note ExCorr here is the photodiode signal.
        elif self.FileType.value > 1:
            self.Trace = [0]*self.NumSamples
        
        self.ReadSpecData()
        return

    def ReadHeaderInfo(self):
        '''
        Read the header (first 7 lines) and extract useful info.
        This is called in initialization.
        
        Useful info that is set:
        - Start date and time (as a time struct).
        - Number of data points acquired.
        - Excitation Wavelength range.
        - Emission Wavelength Range.
        - Run type (from RunTypes enum).
        '''
        #Read the header info to determine the run type
        with file(self.FilePath, 'r') as thefile:
            if self.FileType == self.FileTypes.Session:
                success = self._ReadHdrSession(thefile)
            elif self.FileType == self.FileTypes.Trace:
                success = self._ReadHdrTrace(thefile)
            elif self.FileType == self.FileTypes.Group:
                success = self._ReadHdrGroup(thefile)
        return success
        
    def _ReadHdrSession(self, thefile):
        for i, line in enumerate(thefile):
            if i==1:
                wrds = line.split()
                self.AcqStart = time.strptime(wrds[-2] + ' ' + wrds[-1],
                                              '%Y-%m-%d %H:%M:%S')
            elif i==5:
                wrds = line.split()
                self.NumSamples = int(wrds[0])
            elif i==6:
                success = self._ReadWLRangeLine(line)
                break
        return success

    def _ReadHdrTrace(self, thefile):
        for i, line in enumerate(thefile):
            if i==1:
                #No acquisition time in trace files, just use file creation time as an estimate.
                self.AcqStart = time.localtime(os.path.getctime(self.FilePath))
                self.NumSamples = int(line)
            elif i==2:
                success = self._ReadWLRangeLine(line)
                break
        return success

    def _ReadHdrGroup(self, thefile):
        #No acquisition time in trace files, just use file creation time as an estimate.
        success = True
        self.AcqStart = time.localtime(os.path.getctime(self.FilePath))
        self.PMTmode = 'Correction'
        if 'excorr' in self.FilePath:
            self.RunType = self.RunTypes.Excitation
        elif 'emcorr' in self.FilePath:
            self.RunType = self.RunTypes.Emission
        self.NumSamples = -100
        for i, line in enumerate(thefile):
            if i==3:
                self.NumSamples = int(line)
            elif i==6:
                wrds = line.split()
                if self.RunType == self.RunTypes.Excitation:
                    self.ExRange = [float(wrds[0])]
                elif self.RunType == self.RunTypes.Emission:
                    self.EmRange = [float(wrds[0])]
                else:
                    print("ERROR!! Bad correction file (the names need excorr or emcorr).")
                    success = False
                    break
            elif i==6+self.NumSamples-1:
                wrds = line.split()
                if self.RunType == self.RunTypes.Excitation:
                    self.ExRange += [float(wrds[0])]
                elif self.RunType == self.RunTypes.Emission:
                    self.EmRange += [float(wrds[0])]
                break
        return success
    
    def _ReadWLRangeLine(self, line):
        success = True
        if line[0]=='D':
            self.PMTmode = 'Digital'
        elif line[0] == 'A':
            self.PMTmode = 'Analogue'
        else:
            self.PMTmode = 'Unknown'
            success = False
        wrds = line.split()
        self.ExRange = [float(val) for val in wrds[1].split(':')[0].split('-')]
        self.EmRange = [float(val) for val in wrds[1].split(':')[1].split('-')]
        if len(self.ExRange)>1 and len(self.EmRange)>1:
            self.RunType = self.RunTypes.Synchronous
        elif len(self.ExRange)>1:
            self.RunType = self.RunTypes.Excitation
        elif len(self.EmRange)>1:
            self.RunType = self.RunTypes.Emission
        else:
            self.RunType = self.RunTypes.Unknown
            success = False
        return success

    def ReadSpecData(self):
        '''
        Read the data from the file.
        
        If the file is a trace, this will read:
        - WL (list of wavelengths)
        - Trace (the trace - the caller is expected to know what it is)
        If the file is a session, this will also read:
        - Spec (the spectrum)
        - SpecRaw (uncorrected for excitation/emission)
        - ExCorr (the excitation correction data from the photodiode)
        '''
        if self.FileType == self.FileTypes.Session:
            self._ReadSessionData()
        elif self.FileType == self.FileTypes.Trace:
            self._ReadTraceData()
        elif self.FileType == self.FileTypes.Group:
            self._ReadGroupData()
        return
        
    def _ReadSessionData(self):
        with file(self.FilePath, 'r') as thefile:
            for i, line in enumerate(thefile):
                if i > 7 and i < (8 + self.NumSamples):
                    wrds = line.split()
                    self.WL[i-8] = float(wrds[0])
                    self.SpecRaw[i-8] = float(wrds[1])
                    self.Spec[i-8] = float(wrds[3])
                elif i > (8 + self.NumSamples + 7) and \
                   i < (8 + self.NumSamples + 7 + self.NumSamples):
                    self.ExCorr[i-(8+self.NumSamples+7)] = float(wrds[1])
        return

    def _ReadTraceData(self):
        with file(self.FilePath, 'r') as thefile:
            for i, line in enumerate(thefile):
                if i > 3 and i < (4 + self.NumSamples):
                    wrds = line.split()
                    self.WL[i-4] = float(wrds[0])
                    self.Trace[i-4] = float(wrds[1])
        return

    def _ReadGroupData(self):
        with file(self.FilePath, 'r') as thefile:
            for i, line in enumerate(thefile):
                if i > 5 and i < (6 + self.NumSamples):
                    wrds = line.split()
                    self.WL[i-6] = float(wrds[0])
                    self.Trace[i-6] = float(wrds[1])
        return
