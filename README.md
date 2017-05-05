# ESPECTRO
ESPECTRO is a Python library to facilitate the analysis of high-resolution spectral data reduced with the OPERA pipeline. For more details about data reduction with OPERA see http://wiki.lna.br/wiki/espectro. 

To start using the ESPECTRO tools, download the following libraries:

* `spectralclass.py`
* `espectrolib.py`

Example:
```python
  from spectralclass import Spectrum
  spc = Spectrum("spectrum.m.fits.gz")
  spc.info()
```
One can find more examples to use the ESPECTRO libraries in the App's available (any file starting with `App_`). The user can use these examples as starting point to develop their own applications.  

For example, the application `App_extract.py` extracts the spectrum within a given wavelength range from the OPERA FITS product.  Execute the following command from a Terminal shell:

`
$ESPECTRO_PATH/App_extract.py --input=spectrum.m.fits.gz --wlrange="650 660" --spectype=norm -tr
`

where the input data file `spectrum.m.fits.gz` is a spectrum product reduced by OPERA. In the example it will print out the nomalized (option `--spectype=norm`) spectrum, i.e., the following three columns:

`
wavelength(nm) flux flux_err
`

where the wavelength range will be between 650 and 660 nm (option `--wlrange="650 660"`). The wavelength is already corrected by heliocentric observer's velocity (option `-r`) and the wavelength is also corrected using telluric lines as reference (option `-t`).
