# ESPECTRO
Python library to perform analysis on high-resolution spectral data reduced with OPERA pipeline. For more details about data reduction with OPERA see http://wiki.lna.br/wiki/espectro). 

To start using the ESPECTRO tools, download the following libraries:

* `spectralclass.py`
* `espectrolib.py`

Example:
```python
  from spectralclass import Spectrum
  spc = Spectrum("spectrum.m.fits.gz")
  spc.info()
```
In order to provide examples on the utilities of ESPECTRO libraries, there are a few App's available, so the user can use them as a starting point to develop their own applications.  

For example, the application `App_extract.py` can be used to extract the spectrum within a given wavelength range. 

`
$ESPECTRO_PATH/App_extract.py --input=N20160912G0053.m.fits.gz --wlrange="650 660" --spectype=norm -tr
`

In the example above it will extract the nomalized spectrum (option `--spectype=norm`), within the wavelength range between 650 and 660 nm (option `--wlrange="650 660"`), where the wavelength is already corrected by the heliocentric velocity (option `-r`) and the wavelength is also corrected using telluric lines as reference (option `-t`).
