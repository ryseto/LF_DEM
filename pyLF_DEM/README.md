# pyLFDEM

pyLFDEM is a Python wrapper for [LF_DEM](https://github.com/rmari/LF_DEM).

## Installation

The wrapping is actually done by [SWIG](http://www.swig.org/) (version 3),
resulting in the files `pyLFDEM.cxx` and `pyLFDEM.py`.
You can either use these two files as provided in this repository,
in which case you do not need to worry about SWIG and can skip the first step below,
or you can generate your own wrapping code
(for instance if you want to use another version of `LF_DEM`).

To generate new `pyLFDEM.cxx` and `pyLFDEM.py` files:

`$ swig -python -c++ Simulation.i`

To compile the actual `pyLFDEM`, you need to provide your own environment
variables in a file called `setup_in.py`.
A generic file `setup_in_Linux.py` is provided in the repository.
You can modify this file for your own case and copy it as `setup_in.py`.
Once you have a `setup_in.py`, you just need to run

`$ python setup.py build_ext --inplace`

If successful, the installation creates a shared library
`_pyLFDEM*.so` (the extension is probably OS dependent, it's `.so` on Linux)
and a Python wrapping code `pyLFDEM.py`.
You can copy these two files to a folder within your Python path.

There is a `Makefile` to take care of these 3 steps.
Just provide a file called `Makefile.in` defining your `swig` binary in `SWIG`
and your desired installation directory in `PYTHON_SITE_PACKAGES`.
The full install is then:

`$ make swig; make; make install`


## Usage

To instantiate a Simulation() class from python, just do
```
import pyLFDEM as lf

simu = lf.Simulation()
```

From there you can access all public
methods and members of the Simulation class.
There are not that many but you can
do a lot with it. In particular, the Simulation object exposes its
System and ParameterSet instances!

There are example scripts in the `demo` folder, but these demos
may be outdated. Sorry...
