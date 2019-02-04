# Automated Fragmentation AIMD Calculation

[![python version](https://img.shields.io/pypi/pyversions/aimdfragmentation.svg?logo=python&logoColor=white)](https://pypi.org/project/aimdfragmentation)
[![PyPI](https://img.shields.io/pypi/v/aimdfragmentation.svg)](https://pypi.org/project/aimdfragmentation)
[![Build Status](https://travis-ci.com/njzjz/aimdfragmentation.svg?branch=master)](https://travis-ci.com/njzjz/aimdfragmentation)
[![Coverage Status](https://coveralls.io/repos/github/njzjz/aimdfragmentation/badge.svg?branch=master)](https://coveralls.io/github/njzjz/aimdfragmentation?branch=master)
[![codecov](https://codecov.io/gh/njzjz/aimdfragmentation/branch/master/graph/badge.svg)](https://codecov.io/gh/njzjz/aimdfragmentation)

A automated fragmentation method for Ab Initio Molecular Dynamics (AIMD).

**Author**: Jinzhe Zeng

**Email**: jzzeng@stu.ecnu.edu.cn

## Requirements
* [OpenBabel](https://github.com/openbabel/openbabel)
* [numpy](https://github.com/numpy/numpy)
* [ASE](https://gitlab.com/ase/ase)
* [GaussianRunner](https://github.com/njzjz/GaussianRunner)

## Installation

### Using pip

```sh
$ pip install aimdfragmentation

```

### Build from source

You should install [Gaussian 16](http://gaussian.com/gaussian16/) and [OpenBabel](http://openbabel.org) first. Then:

```sh
git clone https://github.com/njzjz/aimdfragmentation
cd aimdfragmentation/
pip install .
```

## Example

### Run a Python program

You can see [examples/example.py](examples/example.py) as an example, and run with:

```sh
python example.py
```

### Run MD with LAMMPS

See [njzjz/Pyforce](https://github.com/njzjz/Pyforce) repository and install Pyforce module. Then rename [examples/example.py](examples/example.py) as `force.py` and put it where you run LAMMPS. Add a line in the LAMMPS input file:

```
fix 1 all pyforce C H O
```

