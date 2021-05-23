# Automated Fragmentation AIMD Calculation

[![DOI:10.3390/molecules26113120](https://img.shields.io/badge/DOI-10.3390%2Fmolecules26113120-blue.svg)](https://doi.org/10.3390/molecules26113120)
[![python version](https://img.shields.io/pypi/pyversions/aimdfragmentation.svg?logo=python&logoColor=white)](https://pypi.org/project/aimdfragmentation)
[![PyPI](https://img.shields.io/pypi/v/aimdfragmentation.svg)](https://pypi.org/project/aimdfragmentation)
[![codecov](https://codecov.io/gh/njzjz/aimdfragmentation/branch/master/graph/badge.svg)](https://codecov.io/gh/njzjz/aimdfragmentation)

A automated fragmentation method for Ab Initio Molecular Dynamics (AIMD).

Fragment-based Ab Initio Molecular Dynamics Simulation for Combustion, Molecules, 2021, 26 (11), 3120, DOI: [10.3390/molecules26113120](https://doi.org/10.3390/molecules26113120).

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

