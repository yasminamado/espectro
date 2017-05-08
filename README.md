# ESPECTRO
ESPECTRO is a Python library to facilitate the analysis of high-resolution spectral data reduced with the OPERA pipeline. For more details about data reduction with OPERA see http://wiki.lna.br/wiki/espectro. 

To start using the ESPECTRO tools, download the following libraries:

* `spectralclass.py`
* `espectrolib.py`

Usage example:
```python
  from spectralclass import Spectrum
  spc = Spectrum("spectrum.m.fits.gz")
  spc.info()
```
One can find more examples on how to use the ESPECTRO libraries in the App's available (any file starting with `App_`). The user can also use these examples as starting point to develop their own applications.  

As an example, the application `App_extract.py` extracts the spectrum within a given wavelength range from the OPERA FITS product (`*.m.fits.gz`).  Run the following command from a Terminal shell:

`
python $ESPECTRO_PATH/App_extract.py --input=spectrum.m.fits.gz --wlrange="650 660" --spectype=norm -tr -e
`

The input data file `spectrum.m.fits.gz` is a spectrum product reduced by OPERA. In the example above it will print out the nomalized spectrum (option `--spectype=norm`). It will print out the following three columns:

`
wavelength(nm) flux flux_err
`

The wavelength range is selected by the option `--wlrange="650 660"`. The option `-r` applies the heliocentric velocity correction and the option `-t` applies wavelength correction using telluric lines as reference. One may these options to get the spectrum without these correction. Option -e prints an informative header about the spectrum. 
