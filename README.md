# Entropic Trajectories

[![Join the chat at https://gitter.im/EntropicTrajectories/entropictrajectories](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/EntropicTrajectories/entropictrajectories)

This is the [Entropic Trajectories library](https://github.com/EntropicTrajectories/entropictrajectories), which provides a direct wrapping of functionality.  

## Prerequisites
The Entropic Trajectories Framework uses many existing packages under the hood.  These include the BLAS and LAPACK libraries, which must be installed before cloning this repository.  BLAS and LAPACK can be found at,
 - BLAS (Basic Linear Algebra Subprograms); http://www.netlib.org/blas/
 - LAPACK (Linear Algebra PACKage); http://www.netlib.org/lapack/
### Linux installation
To install BLAS and LAPACK on Ubuntu simply run the following command in the terminal,
```
$ sudo apt-get install libblas-dev liblapack-dev
```
Other distributions can likely issue a similar command with different package managers (See https://distrowatch.com/dwres.php?resource=package-management for a list of common ones).  You will also need a fortran compiler if you don't already have one,
```
$ sudo apt-get install gfortran
```
Again, other distributions will use a similar command with a different package manager.

### Windows


## Pybind11
The ET framework uses a python wrapper called *pybind11* in order to generate python bindings.
 - Pybind11; https://github.com/pybind/pybind11
 (Wenzel Jakob)
Pybind is compiled together with the ET framework, so one only needs to extract the source files into the include directory.  This is done automatically via a submodule.


## Installing from PyPI

For 64-bit Linux or Mac systems, instally 'etraj' should just require running:

```
pip install etraj
```

You can then test that it works by running the example above.

## Installing from source

Requirements: You must have CMake>=2.8.12 and a C++11 compatible compiler (GCC>=4.8) to build.

First, you must check out this repository using the recursive command in order to download the latest versions of the included libraries as well,
```
$ git clone --recurse-submodules git@github.com:infophysics/EntropicTrajectories.git
```
If the recursive command failed, you can download the main files and simply update the submodules using,
```
$ git clone git@github.com:infophysics/EntropicTrajectories.git
$ git submodule init
$ git submodule update
```
Then, once all of the submodules are there, simply run the installer,

```
python setup.py install
```

## Usage

Python bindings to the Entropic Trajectories library:

```
import etraj


```

For more examples on possible calls, please see the tests folder.

### Support

* Bugs: Please report bugs to the [issue tracker on Github](https://github.com/) such that we can keep track of them and eventually fix them.  Please explain how to reproduce the issue (including code) and which system you are running on.
* Help: Help can be provided also via the issue tracker by tagging your issue with 'question'
* Contributing:  Please fork this repository then make a pull request.  In this pull request, explain the details of your change and include tests.

## Technical implementation

This package is a [pybind11](https://pybind11.readthedocs.io/en/stable/intro.html) wrapper of [NEST](https://github.com/NESTCollaboration/nest) that uses [TravisCI](https://travis-ci.org) to build binaries using the [manylinux](https://github.com/pypa/python-manylinux-demo) [Docker image](https://www.docker.com).

* Implementation also based on [this](http://www.benjack.io/2018/02/02/python-cpp-revisited.html)

See AUTHORS.md for information on the developers.

## Citation

When you use `nestpy`, please say so in your slides or publications (for publications, see Zenodo link above).  You can mention this in addition to how you cite NEST.  This is important for us being able to get funding to support this project.
