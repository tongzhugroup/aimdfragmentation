# Automated Fragmentation AIMD Calculation

[![DOI:10.26434/chemrxiv.11462160](https://zenodo.org/badge/DOI/10.1039/C9CP05091D.svg)](https://doi.org/10.26434/chemrxiv.11462160)
[![python version](https://img.shields.io/pypi/pyversions/aimdfragmentation.svg?logo=python&logoColor=white)](https://pypi.org/project/aimdfragmentation)
[![PyPI](https://img.shields.io/pypi/v/aimdfragmentation.svg)](https://pypi.org/project/aimdfragmentation)
[![codecov](https://codecov.io/gh/njzjz/aimdfragmentation/branch/master/graph/badge.svg)](https://codecov.io/gh/njzjz/aimdfragmentation)

A automated fragmentation method for Ab Initio Molecular Dynamics (AIMD).

Combustion Driven by Fragment-based Ab Initio Molecular Dynamics Simulation, DOI: 10.26434/chemrxiv.11462160

**Author**: Jinzhe Zeng

**Email**: jzzeng@stu.ecnu.edu.cn

## Installation

### Using pip

```sh
pip install aimdfragmentation

```

You also need to install [Gaussian 16](http://gaussian.com/gaussian16/).

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

