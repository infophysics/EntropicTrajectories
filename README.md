# Entropic Trajectories

[![Join the chat at https://gitter.im/EntropicTrajectories/entropictrajectories] [![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/gitterHQ/gitter)

This is the [Entropic Trajectories library](https://github.com/NESTCollaboration/nest), which provides a direct wrapping of functionality.  

## Installing from PyPI

For 64-bit Linux or Mac systems, instally 'etraj' should just require running:

```
pip install etraj
```

You can then test that it works by running the example above.

## Installing from source

Requirements: You must have CMake>=2.8.12 and a C++11 compatible compiler (GCC>=4.8) to build.

First, you must check out this repository then simply run the installer:

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
