# pyMSHRF - HRF Calculator for GC/MS annotation

# Theoritical Background
`High-resolution filtering (HRF)` approach is both feasible and highly specific toward correct identifications, which bridges the gap between unit resolution GC/MS spectra and accurate mass data dubious. [(ref)](https://doi.org/10.1021/acs.analchem.5b01503)

The isotope pattern calculator module is based on pyISOPACh package. [(ref)](https://github.com/AberystwythSystemsBiology/pyISOPACh)


# Installation

pyMSHRF requires Python 3+ and is unfortunately not compatible with Python 2. If you are still using Python 2, a clever workaround is to install Python 3 and use that instead.

The easiest way of installing pyMSHRF is using ```pip```:

```
pip install pyMSHRF
```
# Usage
## Usage of HRF calculation function
```python
import pyMSHRF
import numpy as np

formula = 'C5H7N3O2'
peaks_query = np.array([[191.09071, 14670.0], [124.05742, 3543.0], [141.09334, 6191.0]], dtype = np.float32)

# Add TMS-derived group if desired
formula_derived = pyMSHRF.derivatization(formula, num_tms=1, num_meox=0)

# Calculate HRF score
HRF_score = pyMSHRF.HRF(formula, peaks_query, delta_da = 0.02)
print(f"HRF score is {HRF_score}.")
```
Calculate the `reverse high-resolution filtering (RHRF)` score need reference spectrum:

```python
peaks_reference = np.array([[82, 6.99], [141, 999.0], [124, 562.49]], dtype = np.float32)

# Calculate RHRF score
RHRF_score = pyMSHRF.RHRF(formula, peaks_query, peaks_reference, delta_da = 0.02)
print(f"RHRF score is {RHRF_score}.")
```


## Useful functions
### Reading Spectra from a File
For ease of use, a function named `read_one_spectrum` is provided in the pyMSHRF package, allowing you to easily read spectra from a file.[(ref)](https://github.com/YuanyueLi/MSEntropy) Here is an example of how you can use it:

```python
import pyMSHRF

# Load all spectra from file into python
spectra_list = list(pyMSHRF.read_one_spectrum('path/to/spectrum/file'))

# Get spectrum peak list
for spectrum in spectra_list:
    query_spec = spectrum['peaks']
```
This function returns a dictionary, where each key-value pair corresponds to a specific metadata of the spectrum.

Currently, the `read_one_spectrum` function supports the following file formats: `.mgf`, `.msp`, `.mzML`, and `.lbm2`.

# Reference
1. N. W. Kwiecien et al., High-Resolution Filtering for Improved Small Molecule Identification via GC/MS. Analytical Chemistry 87, 8328-8335 (2015).

2. <https://github.com/AberystwythSystemsBiology/pyISOPACh>

3. <https://github.com/YuanyueLi/MSEntropy>
