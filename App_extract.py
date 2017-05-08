# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Estract spectrum from OPERA fits products
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    ./extract.py --input=1830317o.pol.fits.gz --wlrange="650 665" --spectype=norm -tr
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

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string', default="")
parser.add_option("-w", "--wlrange", dest="wlrange", help="Output wavelength range (nm)",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-e", action="store_true", dest="header", help="print header", default=False)
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with extract.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output wavelength range: ', options.wlrange
    print 'Spectrum type: ', options.spectype
    print 'Print header: ', options.header
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio


spc = Spectrum(options.input, options.spectype, options.polar, options.telluric, options.helio)

if options.header :
    spc.info()
    print "# wavelength(nm) flux flux_err"

if options.wlrange :
    wl0 = float(options.wlrange.split()[0])
    wlf = float(options.wlrange.split()[1])
else :
    wl0 = spc.wl[0]
    wlf = spc.wl[-1]

wl,flux,fluxerr = spc.extractChunk(wl0,wlf)
chunk = SpectrumChunk(wl,flux,fluxerr)
chunk.printdata()
