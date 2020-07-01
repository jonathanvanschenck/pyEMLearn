# pyEMLearn
A computational E&amp;M library written in python, intended for learning the fundamentals

# Installation
`pyEMLearn` is available through `pip`:
```bash
 $ pip3 install pyEMLearn
```
(or just `pip` if you are on windows, or have aliased your `python3` installation)

# Usage
`pyEMLearn` works by building up a stack of layers, then simulation an EM wave
scattering off the system. You can specify the angle of incidence, as well as
the wavelength of the EM wave and calculate the reflection, transmission and
internal electric field. An example of that is shown in `examples/sparse_reflectance.py`

This package assumes that all length measurement are given in microns. You can
ignore this default, but you will not be able to use built in catalog of materials.

## Materials and Layers
Each `Layer` object in the stack is comprised of a `Materials` object and a thickness
(possibly 0 or infinity). `pyEMLearn` provides a small catalog of materials in the
`pyEMLearn.catalog` submodule, which you can check out, as well as some helper
classes for building your own custom material, using `pyEMLearn.materials`. For
example, if you want to model one of you layers using a Sellmeier model for the
index of fraction, use: `pyEMLearn.materials.SellmeierIndex`.

The only `Layer` object currently intended for users is the base `pyEMLearn.layers.Layer`
class. A list of these are intended to be passed to the system solver object,
but nothing really else is user-intended.

## System Solver
The many object that users will interact with is the `pyEMLearn.layers.System`
class which handles all the simulations and calculation. It takes as its arguments
a list of `Layer` objects, and then exposes a few methods for calculating various
simulation parameters. A common workflow is:

```python
import numpy as np
import pyEMLearn.layers as lay
import pyEMLearn.materials as mat
import pyEMLearn.catalog as cat

# Define the system
system = lay.System([
  lay.Layer(cat.dielectrics.BK7,0.1),      # 100nm of BK7 Glass
  lay.Layer(cat.metals.Ag,0.01),           # 10nm of Silver
  lay.Layer(mat.ConstantIndex(2.1),0.05)   # 50nm of custom material with n=2.1
])

system.compile() # Compile the system to prepare for solving

# Define the solution parameters
wls = np.arange(0.3,0.9,0.01) # pick wavelengths in 10nm steps from 300 to 900nm
aoi = len(wls)*[0.0]          # normal incidence (expects radians)
system.solve(wls,aoi)         # solve the system

# Get the results
R,T = system.get_RT("s")      # get the reflectance and transmission for s-polarized
                              #  light, results are arrays of the same length wls
```

# Development
To install toward developing this package: fork and download this repo from github, open a terminal in the head directory and run:
```bash
 $ python3 -m venv venv
 $ source venv/bin/activate
 (venv) $ pip install -r requirements.txt
 (venv) $ python setup.py install
```
This will create a virtual environment, activate it, install the required packages,
then install the current version of `pyEMLearn` into your virtual environment.

After you have made changes to the source code, re-run the last line to install
your modified version.

## Creating source packages
```bash
 $ python3 setup.py sdist bdist_wheel
```

## Uploading to PyPI
```bash
 $ python3 -m twine upload dist/*
```

## Development Philosophy
This package is really intended as a learning tool, rather than a full computational workhorse. As a result, much of the core code is *deliberately* inefficient so as to be more closely connected with the theory of the scattering (transfer) matrix methods. There are plenty of other computational E&amp;M libraries which are both more efficient and have more sophisticated features (like [`EMpy`](https://github.com/lbolla/EMpy), [`EMUstack`]( https://github.com/bjornsturmberg/EMUstack) and [`Rigorous-Couple-Wave-Analysis`](https://github.com/zhaonat/Rigorous-Coupled-Wave-Analysis)
] to name a few). So, development in `pyEMLearn` will always prioritize user-clarity over code-efficiency.

# Acknowledgements
The authors would like to thank Dr. Raymond Rumpf (of UTEP) for his extensive [educational resources](https://empossible.net/academics/emp5337/), on which this project is principally based.

# To Do
 - [ ] Comment Everything . . .
 - [ ] Add some descriptive examples
 - [ ] Add some better examples
 - [ ] Restructure the catalog
 - [ ] Add more materials to the catalog
 - [x] Add transmittance matrix approach
 - [x] Create `source` object, which can produce k-vecs and polarization vecs
 - [x] Create `field` object? which can track E?
 - [ ] Add Poynting vector functionality to `field` object
 - [x] Begin thinking about RCWA
 - [ ] Fix RCWA, its broke AF
