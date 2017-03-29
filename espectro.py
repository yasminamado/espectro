#!/home/yasmin/anaconda2/bin/python
# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
            #!/Users/edermartioli/Local/Ureka/variants/common/bin/python
    Created on Mar 10 2017
    
    Description: A wrapper to run ESPECTRO V1.0
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    ./espectro.py --input=1830317o.pol.fits.gz --wlrange="650 665"
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
import espectrolib

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string',default="")
parser.add_option("-w", "--wlrange", dest="wlrange", help="Output wavelength range (nm)",type='string',default="")
parser.add_option("-p", action="store_true", dest="polar", help="polar",default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with espectro.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output wavelength range: ', options.wlrange
    print 'Polar? ', options.polar

wl0, wlf = options.wlrange.split()

spc = Spectrum(options.input, polar=options.polar)

#spc.info()
wl,flux = spc.extractChunk(float(wl0), float(wlf))

dl = 2.0
wlranges = [(float(wl0),float(wl0)+dl),(float(wlf)-dl,float(wlf))]
wl,flux = espectrolib.removeBackground(wl,flux,wlranges)

for i in range(len(wl)) :
    print wl[i], flux[i]

#spc.printdata()
