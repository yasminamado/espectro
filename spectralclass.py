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
from scipy.interpolate import UnivariateSpline
from scipy import stats

########## SPECTRUM CLASS ############
class Spectrum :
    'Common base class for a spectrum'

    def __init__(self, Filename, Spectype="raw", polar=False, Telluric=False, Helio=False, SourceRV=0.0):

        """
        Create a Spectrum object.
        
        Parameters
        ----------
        filename : string
            File to read the spectrum from.

        Examples
        --------
        >>> spc = Spectrum("spectrum.m.fits.gz")
        """
        self.sourceRV = SourceRV
        
        self.telluric = Telluric
        self.helio = Helio
        self.spectype = Spectype
        self.filepath = Filename
        basename = espectrolib.getbasename(self.filepath)
        
        self.id = basename[0]
        self.filename = os.path.basename(self.filepath)

        try :
            if self.filepath.endswith(".m.fits.gz") or self.filepath.endswith(".m.fits"):
                self.wl,self.flux,self.fluxerr=self.loadSpectrumFromMFITS(self.filepath, self.spectype, self.telluric, self.helio)
            elif self.filepath.endswith(".pol.fits.gz") or self.filepath.endswith(".pol.fits") :
                self.wl,self.flux,self.fluxerr=self.loadSpectrumFromPOLFITS(self.filepath, self.spectype, polar, self.telluric, self.helio)
            else :
                print "Error: file type not supported for input spectrum: ",self.filepath
                exit()
        except :
            print "Error: could not open file: ",self.filepath
            exit()

    #--- Function to load spectrum from .m.fits.gz file
    def loadSpectrumFromMFITS(self,fitsfilename, spectype="raw", telluric=False, helio=False, sort=True):

        wl,flux = [],[]
        hdu = fits.open(fitsfilename)
        
        self.header = hdu[0].header
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']
        self.HJDTT = hdu[0].header['HJDTT']
        self.exptime = hdu[0].header['EXPTIME']

        self.orders = hdu[0].data[0]

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

        if sort :
            indices = wltmp.argsort()
            wl = wltmp[indices]
            flux = influx[indices]
            fluxerr = influxerr[indices]
            self.orders = self.orders[indices]
        else :
            wl = wltmp
            flux = influx
            fluxerr = influxerr

        return wl,flux,fluxerr
    #------------

    #--- Function to load spectrum from .pol.fits.gz file
    def loadSpectrumFromPOLFITS(self,fitsfilename, spectype="raw", polar=False, telluric=False, helio=False, sort=True):
    
        wl,flux = [],[]
        hdu = fits.open(fitsfilename)
        self.header = hdu[0].header
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']
        self.HJDTT = hdu[0].header['HJDTT']
        self.exptime = hdu[0].header['EXPTIME']

        self.orders = hdu[0].data[0]

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
        
        if sort :
            indices = wltmp.argsort()
            wl = wltmp[indices]
            flux = influx[indices]
            fluxerr = influxerr[indices]
            self.orders = self.orders[indices]
        else :
            wl = wltmp
            flux = influx
            fluxerr = influxerr

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
        print "#**************************"
        print "#Info for spectrum: ",self.filename
        print "#Object:", self.object
        print "#Instrument:",self.instrument
        print "#wl0,wlf =",self.wl[0],",",self.wl[-1],"nm"
        sampling = (self.wl[-1] - self.wl[0])/float(len(self.wl))
        print "#sampling =",sampling," nm/pixel"
        print "#Wave telluric correction =",self.telluric
        print "#Wave heliocentric correction =",self.helio
        print "#Flux type =",self.spectype
        print "#HJD (TT) =",self.HJDTT
        print "#EXPTIME = ",self.exptime, "s"
        print "#<F> =",self.flux.mean(),"+-",self.flux.std()
        print "#**************************"
    #------------

    #--- Print spectrum data
    def printdata(self) :
        for i in range(len(self.wl)) :
            print self.wl[i],self.flux[i],self.fluxerr[i]
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

    #--- Mask data
    def maskdata(self, lines, width) :
        for line in lines :
            wl0 = line - width
            wlf = line + width
            mask = np.where(np.logical_or(self.wl < wl0, self.wl > wlf))
            self.wl = self.wl[mask]
            self.flux = self.flux[mask]
            self.fluxerr = self.fluxerr[mask]
    #------------

    #--- bin spectrum
    def binning(self, rvsampling_kps, wl0=0.0, wlf=0.0, median=False) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]

        fluxvar = self.fluxerr*self.fluxerr
        
        bins = []
        wl = wl0
        while wl < wlf :
            bins.append(wl)
            wl *= (1.0 + rvsampling_kps*1000/constants.c)
        bins = np.array(bins)
        
        digitized = np.digitize(self.wl, bins)

        wl_new = []
        flux_new = []
        fluxvar_new = []

        for i in range(1, len(bins)):
            if len(self.wl[digitized == i]) :
                wl_new.append(self.wl[digitized == i].mean())
                if median :
                    flux_new.append(np.median(self.flux[digitized == i]))
                else :
                    flux_new.append(self.flux[digitized == i].mean())
                fluxvar_new.append(self.flux[digitized == i].std())

        self.wl = np.array(wl_new)
        self.flux = np.array(flux_new)
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------
    
    #--- Extract order
    def extractOrder(self, order, printdata=False) :
        ordermask = np.where(self.orders == order)
        
        return self.wl[ordermask],self.flux[ordermask],self.fluxerr[ordermask]
    #--------------------------

    #--- filter by orders
    def getMinMaxOders(self) :
        
        minorder = np.min(self.orders)
        maxorder = np.max(self.orders)

        return minorder, maxorder
    #--------------------------

    #--- apply mask (keep only data within the mask)
    def applyMask(self) :
        
        mask = {28: [791.40,820.30],
                29: [765.70,791.50],
                30: [740.60,765.80],
                31: [717.50,741.00],
                32: [695.50,717.80],
                33: [675.00,696.00],
                34: [655.80,675.60],
                35: [637.00,656.40],
                36: [619.50,638.00],
                37: [603.00,620.50],
                38: [586.50,604.50],
                39: [570.50,588.00],
                40: [556.50,572.70],
                41: [534.40,560.20],
                42: [530.00,546.50],
                43: [517.00,534.30],
                44: [505.50,521.20],
                45: [496.00,510.20],
                46: [484.00,499.40],
                47: [473.50,488.80],
                48: [464.00,478.60],
                49: [454.50,468.80],
                50: [445.50,459.20],
                51: [437.30,449.50],
                52: [429.00,440.80],
                53: [423.00,431.40],
                54: [414.00,424.00],
                55: [404.00,416.00]}
        
        for item in mask.items() :
            o = item[0]
            wl0 = mask[item[0]][0]
            wlf = mask[item[0]][1]
            logi_wl = np.logical_or(self.wl < wl0, self.wl > wlf)
            logi_order =  np.logical_and(self.orders>o-1,self.orders<o+1)
            cutmask = np.where(np.logical_not(np.logical_and(logi_wl,logi_order)))
            self.wl = self.wl[cutmask]
            self.flux = self.flux[cutmask]
            self.fluxerr = self.fluxerr[cutmask]
            self.orders = self.orders[cutmask]
    #--------------------------

    #--- save spectrum to output file
    def saveToFile(self, output, format='fits') :

        if format == 'fits' :
            self.header['COMMENT'] = "Data processed by ESPECTRO 1.0"
            self.header['ORIGDATA'] = self.id
            self.header['TELLCORR'] = self.telluric
            self.header['HERVCORR'] = self.helio
            self.header['SPECTYPE'] = self.spectype
            self.header['SOURCERV'] = self.sourceRV

            self.header.set('COL1','Wavelength', 'wavelength (nm)')
            self.header.set('COL2','Flux', 'flux (electron)')
            self.header.set('COL3','FluxError', 'flux error (electron)')

            self.header.remove('COL4')
            self.header.remove('COL5')
            self.header.remove('COL6')
            self.header.remove('COL7')
            self.header.remove('COL8')
            self.header.remove('COL9')
            self.header.remove('COL10')
            self.header.remove('COL11')
            self.header.remove('COL12')
            self.header.remove('COL13')
            self.header.remove('COL14')
            
            scidata = []
            scidata.append(self.wl)
            scidata.append(self.flux)
            scidata.append(self.fluxerr)
            
            if os.path.exists(output) :
                os.remove(output)
            fits.writeto(output, scidata, self.header)
        else :
            print 'file format ',format, ' not supported'
            exit()
    #--------------------------

    #--- calculate continuum
    def continuum(self, binsize = 200, overlap = 100, sigmaclip = 3.0, window = 3) :

        nbins = len(self.wl)/binsize

        wlbin = []
        fluxbin = []

        for i in range(nbins):
    
            idx0 = i*binsize - overlap
            idxf = (i+1)*binsize + overlap
    
            if idx0 < 0 : idx0=0
            if idxf > len(self.wl) : idxf=len(self.wl)-1
    
            wltmp = self.wl[idx0:idxf]
            fluxtmp = self.flux[idx0:idxf]
    
            wlmean = np.mean(wltmp)
    
            wlbin.append(wlmean)
            medflux = np.median(fluxtmp)
            medfluxdev = np.median(np.abs(fluxtmp - medflux))
            filtermask = np.where(np.logical_and(fluxtmp > medflux, fluxtmp < medflux + sigmaclip*medfluxdev))
    
            fluxbin.append(np.max(fluxtmp[filtermask]))

        newwlbin = []
        newfluxbin = []

        newwlbin.append(self.wl[0])
        newfluxbin.append(fluxbin[0])

        for i in range(nbins):
            idx0 = i - window
            idxf = i + 1 + window
    
            if idx0 < 0 : idx0=0
            if idxf > nbins : idxf=nbins-1
    
            slope, intercept, r_value, p_value, std_err = stats.linregress(wlbin[idx0:idxf], fluxbin[idx0:idxf])
    
            newwlbin.append(wlbin[i])
            newfluxbin.append(intercept + slope*wlbin[i])


        s = UnivariateSpline(newwlbin, newfluxbin, s=1)
        self.continuum = s(self.wl)
    #--------------------------

    #--- normalize by continuum
    def normalizeByContinuum(self) :
        try :
            self.flux = self.flux / self.continuum
        except :
            print "Could not perfom continuum normalization."
    #--------------------------

    #--- save continuum spectrum to output file
    def saveContinuumToFile(self, output, format='fits') :
    
        if format == 'fits' :
            self.header['COMMENT'] = "Data processed by ESPECTRO 1.0"
            self.header['ORIGDATA'] = self.id
            
            self.header.set('COL1','Wavelength', 'wavelength (nm)')
            self.header.set('COL2','Flux', 'flux (electron)')
            self.header.set('COL3','FluxError', 'flux error (electron)')
            
            self.header.remove('COL4')
            self.header.remove('COL5')
            self.header.remove('COL6')
            self.header.remove('COL7')
            self.header.remove('COL8')
            self.header.remove('COL9')
            self.header.remove('COL10')
            self.header.remove('COL11')
            self.header.remove('COL12')
            self.header.remove('COL13')
            self.header.remove('COL14')
            
            scidata = []
            scidata.append(self.wl)
            scidata.append(self.continuum)
            scidata.append(self.fluxerr)
            
            if os.path.exists(output) :
                os.remove(output)
            fits.writeto(output, scidata, self.header)
        else :
            print 'file format ',format, ' not supported'
            exit()
    #--------------------------

########## SPECTRUM CHUNK CLASS ############
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

    #----------- get chunk spectrum -----------
    def getSpectrum(self) :
        return self.wl, self.flux, self.fluxerr
    #-----------

    #----------- print data and model -----------
    def printdataWithModel(self) :
        for i in range(len(self.wl)) :
            print self.wl[i], self.flux[i], self.fluxerr[i], self.model[i]
    #-----------

    #----------- print data -----------
    def printdata(self) :
        for i in range(len(self.wl)) :
            print self.wl[i], self.flux[i], self.fluxerr[i]
    #-----------

    #--- save spectrum to output file
    def saveToFile(self, output, format='fits') :
        if format == 'fits' :
            scidata = []
            scidata.append(self.wl)
            scidata.append(self.flux)
            scidata.append(self.fluxerr)

            hdu = fits.PrimaryHDU(scidata)
            hdulist = fits.HDUList([hdu])
            
            hdulist[0].header.set('COL1','Wavelength', 'wavelength (nm)')
            hdulist[0].header.set('COL2','Flux', 'flux (electron)')
            hdulist[0].header.set('COL3','FluxError', 'flux error (electron)')
            hdulist[0].header['COMMENT'] = "Data processed by ESPECTRO 1.0"
            if os.path.exists(output) :
                os.remove(output)
            hdulist.writeto(output)
        else :
            print 'file format ',format, ' not supported'
            exit()
    #--------------------------
