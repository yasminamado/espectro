# ESPECTRO
Python library to perform analysis on high-resolution spectral data reduced with OPERA pipeline (http://wiki.lna.br/wiki/espectro). 

Example:

from spectralclass import Spectrum
spc = Spectrum("spectrum.m.fits.gz")
spc.info()
