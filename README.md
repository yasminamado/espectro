# ESPECTRO
Python library to perform analysis on high-resolution spectral data reduced with OPERA pipeline. For more details about data reduction with OPERA see http://wiki.lna.br/wiki/espectro). 

To start using the ESPECTRO tools, download the following libraries:

* spectralclass.py
* espectrolib.py

Example:

```python
  from spectralclass import Spectrum
  spc = Spectrum("spectrum.m.fits.gz")
  spc.info()
```
