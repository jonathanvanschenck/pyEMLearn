# pyEMLearn
A computational E&amp;M library written in python, intended for learning the fundamentals

# Installation
`pyEMLearn` is available through `pip`:
```bash
 $ pip3 install pyEMLearn
```
(or just `pip` if you are on windows, or have aliased your `python3` installation)

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
# V1 To Do, Obsolete
 - [x] Flesh out interface object with "t" matrix
 - [x] Add interface object to layers objects (and halfplane objects)
 - [ ] Modify utils objects to be more general
 - [ ] Change materials repr functions to NOT return parameters, put that in a "summary" function
 - [ ] Finish Docstring-ing
 - [x] Restructure the materials catalog
 - [ ] Add more asserts to System.__init__ to make sure that transmission layers arn't actually injection layers, or gap layers arn't actually just thickness zero layers, etc.
 - [x] Improve computational efficiency for LHI materials? (using 2c.pdf)
 - [ ] Spend some time making all the class variable names consistent in their case/-/_ usage
 - [x] Add config and license files
 - [x] Get Setup on PyPi
 - [ ] Add anisotropy functionality
 - [ ] Make some better examples (bragg reflectors)
 - [ ] Add parallelization options for parameter sweep calculations
 - [x] Add ellipsometric variable support
 - [x] Add internal field support
