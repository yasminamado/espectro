# -*- coding: utf-8 -*-

"""
Spectral Classes
---------------------------
Created on Mar 10 2017
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""
import os
import numpy as np
from astropy.io import fits
import espectrolib
from scipy import constants
from scipy.optimize import curve_fit

########## SPECTRUM CLASS ############
class Spectrum :
    'Common base class for a spectrum'

    def __init__(self, Filename, spectype="raw", polar=False, telluric=False, helio=False):

        """
        Create a Spectrum object.
        
        Parameters
        ----------
        filename : string
            File to read the spectrum from.

        Examples
        --------
        >>> spc = Spectrum("spectrumfile.spc.gz")
        >>> spc = Spectrum("spectrumfile.m.fits.gz")
        """
        self.sourceRV = 0.0
        
        self.filepath = Filename
        basename = espectrolib.getbasename(self.filepath)
        
        self.id = basename[0]
        self.filename = os.path.basename(self.filepath)

        #try :
        if self.filepath.endswith(".m.fits.gz") or self.filepath.endswith(".m.fits"):
            self.wl,self.flux,self.fluxerr=self.loadSpectrumFromMFITS(self.filepath, spectype, telluric, helio)
        elif self.filepath.endswith(".pol.fits.gz") or self.filepath.endswith(".pol.fits") :
            self.wl,self.flux,self.fluxerr=self.loadSpectrumFromPOLFITS(self.filepath, spectype, polar, telluric, helio)
        else :
            print "Error: file type not supported for input spectrum: ",self.filepath
            exit()
        """
        except :
        print "Error: could not open file: ",self.filepath
        exit()
        """

    #--- Function to load spectrum from .m.fits.gz file
    def loadSpectrumFromMFITS(self,fitsfilename, spectype="raw", telluric=False, helio=False):
        
        wl,flux = [],[]
        hdu = fits.open(fitsfilename)
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']
    
        if spectype == "raw" :
            influx = hdu[0].data[8]
            influxvar = hdu[0].data[9]
        elif spectype == "norm" :
            influx = hdu[0].data[10]
            influxvar = hdu[0].data[11]
        elif spectype == "fcal" :
            influx = hdu[0].data[12]
            influxvar = hdu[0].data[13]
        else :
	    print "Error: spectrum type not recognized: ",spectype
	    exit()

        influxerr = np.sqrt(influxvar)


        if telluric :
            wltmp = hdu[0].data[5]
        else :
            wltmp = hdu[0].data[4]

        if helio :
            wltmp += hdu[0].data[6]

	wltmp *= (1.0 - self.sourceRV*1000.0/constants.c)

        NAN_mask = np.where(np.logical_not(np.isnan(influx)))
            
        indices = (wltmp[NAN_mask]).argsort()
        wl = wltmp[NAN_mask][indices]
        flux = influx[NAN_mask][indices]
	fluxerr = influxerr[NAN_mask][indices]
        
        return wl,flux,fluxerr
    #------------

    #--- Function to load spectrum from .pol.fits.gz file
    def loadSpectrumFromPOLFITS(self,fitsfilename, spectype="raw", polar=False, telluric=False, helio=False):
    
        wl,flux = [],[]
        hdu = fits.open(fitsfilename)
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']

        if polar :
            if spectype == "raw" or spectype == "fcal":
	        influx = hdu[0].data[13]
                influxvar = hdu[0].data[14]
            elif spectype == "norm":
	        influx = hdu[0].data[15]
                influxvar = hdu[0].data[16]
            else :
	       print "Error: spectrum type not recognized: ",spectype
	       exit()
        else :
            if spectype == "raw" :
                influx = hdu[0].data[7]
                influxvar = hdu[0].data[8]
            elif spectype == "norm" :
                influx = hdu[0].data[9]
                influxvar = hdu[0].data[10]
            elif spectype == "fcal" :
                influx = hdu[0].data[11]
                influxvar = hdu[0].data[12]
            else :
	       print "Error: spectrum type not recognized: ",spectype
	       exit()

        influxerr = np.sqrt(influxvar)

        if telluric :
            wltmp = hdu[0].data[4]
        else :
            wltmp = hdu[0].data[3]

        if helio :
            wltmp += hdu[0].data[5]

        wltmp *= (1.0 - self.sourceRV*1000.0/constants.c)
            
        NAN_mask = np.where(np.logical_not(np.isnan(influx)))
            
        indices = (wltmp[NAN_mask]).argsort()
        wl = wltmp[NAN_mask][indices]
        flux = influx[NAN_mask][indices]
	fluxerr = influxerr[NAN_mask][indices]
        
        return wl, flux, fluxerr
    #------------

    #--- resampling spectrum
    def resampling(self, wlsampling, wl0, wlf) :
        npoints = int((wlf-wl0)/wlsampling)
        wl_new = np.linspace(wl0, wlf, npoints)
        flux_new = np.interp(wl_new, self.wl, self.flux)
        self.wl = wl_new
        self.flux = flux_new
    #------------

    #--- Print spectrum information
    def info(self) :
        print "**************************"
        print "Info for spectrum: ",self.filename, " Object:", self.object
        print "Instrument:",self.instrument
        print "wl0 =",self.wl[0],"A -- wlf =",self.wl[-1],"A"
        sampling = (self.wl[-1] - self.wl[0])/float(len(self.wl))
        print "sampling =",sampling," A/pixel"
        print "<F> =",self.flux.mean(),"+-",self.flux.std()
        print "**************************\n"
    #------------

    #--- Print spectrum data
    def printdata(self) :
        for i in range(len(self.wl)) :
            print self.wl[i],self.flux[i]
    #------------

    #--- Extract spectral range
    def extractChunk(self, wl0, wlf, printdata=False) :
        mask = np.where(np.logical_and(self.wl > wl0, self.wl < wlf))

        if printdata :
            wlout = self.wl[mask]
            fluxout = self.flux[mask]
            fluxerrout = self.fluxerr[mask]
            for i in range(len(wlout)) :
                print wlout[i], fluxout[i], fluxerrout[i]

        return self.wl[mask],self.flux[mask],self.fluxerr[mask]
    #------------

    #--- Get time from img header
    def getTimeHJDTT(self) :
        header = fits.getheader(self.filepath,0)
        return header['HJDTT']
    #------------

########## SPECTRUM CLASS ############
class SpectrumChunk :
    'Common base class for a spectrum chunk'

    def __init__(self, Wl, Flux, Fluxerr):
        self.wl,self.flux,self.fluxerr = Wl, Flux, Fluxerr


    def removeBackground(self, backgroundsize, wl0=0.0, wlf=0.0) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]
        wlranges = [(float(wl0),float(wl0)+backgroundsize),(float(wlf)-backgroundsize,float(wlf))]
        self.wl,self.flux = espectrolib.removeBackground(self.wl,self.flux,wlranges)


    def snrfilter(self,snrcut) :
        snr = self.flux / self.fluxerr 
        mask = np.where(snr > snrcut)
	self.wl,self.flux,self.fluxerr = self.wl[mask],self.flux[mask],self.fluxerr[mask]


    #--- resampling spectrum
    def resampling(self, wlsampling, wl0=0.0, wlf=0.0) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]

        npoints = int((wlf-wl0)/wlsampling)
        wl_new = np.linspace(wl0, wlf, npoints)
        flux_new = np.interp(wl_new, self.wl, self.flux)
        fluxvar = self.fluxerr*self.fluxerr
        fluxvar_new = np.interp(wl_new, self.wl, fluxvar)

        self.wl = wl_new
        self.flux = flux_new
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------


    #--- bin spectrum
    def binning(self, wlsampling, wl0=0.0, wlf=0.0, median=False) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]

        fluxvar = self.fluxerr*self.fluxerr

	npoints = int((wlf-wl0)/wlsampling)
        bins = np.linspace(wl0, wlf, npoints)
	digitized = np.digitize(self.wl, bins)

	wl_new = [self.wl[digitized == i].mean() for i in range(1, len(bins))]

        if median :
            flux_new = [np.median(self.flux[digitized == i]) for i in range(1, len(bins))]
        else :
            flux_new = [self.flux[digitized == i].mean() for i in range(1, len(bins))]

        fluxvar_new = [self.flux[digitized == i].std() for i in range(1, len(bins))]

        self.wl = wl_new
        self.flux = flux_new
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------


    #----------- FFT filtering -----------
    def fft_filter(self, step=0.01, samp=5):
        """
            Perform Fast Fourier Transform filtering of flux data
        
            Parameters
            ----------
            step : step to define frequency sampling
            samp : number of sample points 
        """
        fourier = np.fft.fft(self.flux)
        n = len(fourier)
        freq = np.fft.fftfreq(n, d=step)
        iy = np.zeros(n, dtype=complex)
    
        for j in xrange(n):
            if -samp < freq[j] < samp:
                iy[j] = fourier[j]
            else:
                iy[j] = 0
        self.flux = np.real(np.fft.ifft(iy, n))
    #-----------

    #----------- Fit Gaussian -----------
    def fitgaussian(self, guess) :

        popt, pcov = curve_fit(espectrolib.gaussfunc, self.wl, self.flux, p0=guess)

        self.model = espectrolib.gaussfunc(self.wl, *popt)

        wlcen = []
        amp = []
        sigma = []

        for i in range(0, len(popt), 3) :
            wlcen.append(popt[i])
            amp.append(popt[i+1])
            sigma.append(popt[i+2])
 
        return wlcen, amp, sigma
    #-----------

    def printdataWithModel(self) :
        for i in range(len(self.wl)) :
            print self.wl[i], self.flux[i], self.fluxerr[i], self.model[i]

    def printdata(self) :
        for i in range(len(self.wl)) :
            print self.wl[i], self.flux[i], self.fluxerr[i]

