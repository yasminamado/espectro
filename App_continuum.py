# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Calculate continuum spectrum
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_continuum.py --input=1830317o.pol.fits.gz --spectype=fcal -t -r
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
from spectralclass import SpectrumChunk
import espectrolib
import numpy as np

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_continuum.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

spectrallines = [384,389,397,403,410,427,434,439,448,454,465,469,471,486,492,502,517,541,587,656,668,707,724,728,824,845,846,850,854,860,866,875,886,901,923,955,1012]
telluriclines = [762, 765, 942]

lines = np.concatenate((spectrallines,telluriclines),axis=0)

spc = Spectrum(options.input, options.spectype, options.polar, options.telluric, options.helio)

width = 2.0
spc.maskdata(lines,width)

continuum_wlsampling = 10.0
spc.binning(continuum_wlsampling, median=True)

spc.printdata()




