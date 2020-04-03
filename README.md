# pyEMLearn
A computational E&amp;M library written in python, intended for learning the fundamentals


# Development
To install toward developing this package: fork and download this repo from github, open a terminal in the head directory and run:
```bash
python3 config.py install
```
(note, if you are on windows or using something like conda which aliases python3, substitute `python` for `python3`)

When you have finished developing, create a pull request through github.

## Development Philosophy
This package is really intended as a learning tool, rather than a full computational workhorse. As a result, much of the core code is *deliberately* inefficient so as to be more closely connected with the theory of the scattering (transfer) matrix methods. There are plenty of other computational E&amp;M libraries which are both more efficient and have more sophisticated features (like [EMpy](https://github.com/lbolla/EMpy), [EMUstack]( https://github.com/bjornsturmberg/EMUstack) and [Rigorous-Couple-Wave-Analysis](https://github.com/zhaonat/Rigorous-Coupled-Wave-Analysis)
] to name a few). So, development in pyEMLearn will always prioritize user-clarity over code-efficiency.

# Acknowledgements
The authors would like to thank Dr. Raymond Rumpf (of UTEP) for his extensive [educational resources](https://empossible.net/academics/emp5337/), on which this project is principally based.

# To Do
* Flesh out interface object with "t" matrix
* Add interface object to layers objects (and halfplane objects)
* Modify utils objects to be more general
* Change materials repr functions to NOT return parameters, put that in a "summary" function
* Finish Docstring-ing
* Restructure the materials catalog
* Add more asserts to System.__init__ to make sure that transmission layers arn't actually injection layers, or gap layers arn't actually just thickness zero layers, etc.
* Improve computational efficiency for LHI materials? (using 2c.pdf)
* Spend some time making all the class variable names consistent in their case/-/_ usage
* Add config and liscene files
* Get Setup on PyPi
* Restructure and flesh-out the material database
* Add anisotropy functionality
* Change examples over to jupyter
* Make some better examples (bragg reflectors)
* Think about PWA and RCWA as extensions
* Add parallelization options for parameter sweep calculations
* Add ellipsometric variable support
* Add internal field support