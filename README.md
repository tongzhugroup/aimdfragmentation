# Automated Fragmentation AIMD Calculation
[![python3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://badge.fury.io/py/AIMDFragmentation)[![pypi](https://badge.fury.io/py/AIMDFragmentation.svg)](https://badge.fury.io/py/AIMDFragmentation)

A automated fragmentation method for Ab Initio Molecular Dynamics (AIMD).

**Author**: Jinzhe Zeng
**Email**: jzzeng@stu.ecnu.edu.cn

## Requirements
* [numpy](https://github.com/numpy/numpy)
* [GaussianRunner](https://github.com/njzjz/GaussianRunner)

## Installation
### Using pip
```sh
$ pip install AIMDFragmentation
```
### Build from source
You should install [Gaussian 16](http://gaussian.com/gaussian16/) and [OpenBabel](http://openbabel.org) first. Then build AIMDBlock:
```sh
$ git clone https://github.com/njzjz/AIMDFragmentation.git
$ cd AIMDFragmentation/
$ python setup.py install
```
## Example

### Run a Python program
You can see [examples/example.py](examples/example.py) as an example, and run with:
```sh
python3 example.py
```

### Run MD with LAMMPS
See [njzjz/Pyforce](https://github.com/njzjz/Pyforce) repository and install Pyforce module. Then rename [examples/example.py](examples/example.py) as `force.py` and put it where you run LAMMPS. Add a line in the LAMMPS input file:
```
fix 1 all pyforce C H O
```

