# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 2017
@author: Eder Martioli
Description: ESPECTRO Library
Laboratorio Nacional de Astrofisica, Brazil
"""

import os
import numpy as np
from astropy.io.fits import getheader

######################
def get_fitsfilepaths(directory):
    
    """
    Generates a list of file names in a directory tree
    by walking the tree either top-down or bottom-up.
        
    Parameters
    ----------
    directory : directory path
        
    Returns
    -------
    file_paths: a list of file paths
    """
        
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree:
    for root, directories, files in os.walk(directory):
        for filename in files:
			# Merge strings to form full filepaths
            filepath = os.path.join(root, filename)
            if filename.endswith(".fits") :
                file_paths.append(filepath)

    file_paths.sort()
    return file_paths
######################

############# get basename from FITS file path ###############
def getbasename(filepath) :
    
    basename = ''
    base = os.path.basename(filepath)
    
    dir = os.path.dirname(filepath)
    
    if ".fits.gz" in base :
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    elif ".fits" in base :
        basename = os.path.splitext(base)[0]
    elif ".txt" in base :
        basename = os.path.splitext(base)[0]
    elif ".spc.gz" in base :
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    else :
        print "Error: unknown extension in file: ", filepath
        exit()

    return basename, dir
###########################################

############# get Heliocentric Radial Velocity from header ####
def getEspaonsHelioRV(header) :
    comments = header['COMMENT']
    rv = 0.0
    for c in comments :
        if "Heliocentric velocity" in c :
            rv = float(c[c.find(": ")+1:c.find(" km")])
            break
    return rv
###########################################

############# FFT filtering ###############
def fft_filter(y, step=0.01, samp=20):
    """
        Perform Fast Fourier Transform filtering of data
        
        Parameters
        ----------
        y : data array
        step : step to define frequency sampling
        samp : number of sample points

        Returns
        -------
        yfftclean: filtered data array
    """

    fourier = np.fft.fft(y)
    n = len(fourier)
    freq = np.fft.fftfreq(n, d=step)
    iy = np.zeros(n, dtype=complex)
    
    
    for j in xrange(n):
        if -samp < freq[j] < samp:
            iy[j] = fourier[j]
        else:
            iy[j] = 0
    yfftclean = np.real(np.fft.ifft(iy, n))
    
    return yfftclean
#####

############# Remove Background  ###############
def removeBackground(wl,flux,wlranges):
    """
        Remove background of data
        
        Parameters
        ----------
        wl : wavelength data array
        flux : flux data array
        wlranges : wavelength ranges where to estimate background
        
        Returns
        -------
        wl,flux: reduced data arrays
    """

    mask1 = np.where(np.logical_and(wl > wlranges[0][0], wl < wlranges[0][1]))
    mask2 = np.where(np.logical_and(wl > wlranges[1][0], wl < wlranges[1][0]))
    
    flux_bkg = np.append(flux[mask1],flux[mask2])
    
    background = np.median(flux_bkg)
    
    flux -= background
    
    return wl,flux
#####

############# Multiple gaussian function  ###############
def gaussfunc(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y
#####


######################
def generateList(directory, objectname, objectkey="OBJECT"):
    
    """
    Generates a list of file names in a directory tree
    by walking the tree either top-down or bottom-up.
        
    Parameters
    ----------
    directory : directory path
        
    Returns
    -------
    file_paths: a list of file paths
    """
        
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree:
    for root, directories, files in os.walk(directory):
        for filename in files:
	    # Merge strings to form full filepaths
            filepath = os.path.join(root, filename)
            
            if filename.endswith("m.fits.gz") or filename.endswith("pol.fits.gz") \
               or filename.endswith("m.fits") or filename.endswith("pol.fits"):
                header = getheader(filepath.rstrip('\n'),0)
                objname = header[objectkey]
                if objname == objectname :
                    file_paths.append(filepath)

    file_paths.sort()
    return file_paths
######################


######################
def wlrange(wlstr, spc):
    
    """
        Get wavelength range from input command line string string
        
        Parameters
        ----------
        wlstr : input range in string format: "wl0 wlf". E.g. "500 550"
        spc   : Spectrum() class
        
        Returns
        -------
        wavelength floats: wl0, wlf
        """

    wl0=0.0
    wlf=0.0
    
    if wlstr :
        wl0 = float(wlstr.split()[0])
        wlf = float(wlstr.split()[1])
    else :
        wl0 = spc.wl[0]
        wlf = spc.wl[-1]
    
    return wl0, wlf
######################
