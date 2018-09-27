# Automated Fragmentation AIMD Calculation
A automated fragmentation method for Ab Initio Molecular Dynamics (AIMD).

# Requirements
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
$ python3 setup.py install
```
## Examples
Add the directory of [examples/readbond.py](examples/readbond.py) into PATH and run [readbond.py](examples/readbond.py) directly:
```sh
readbond.py 28 ch4.xyz
```
Forces of [ch4.xyz](examples/ch4.xyz) will be calculated.

Or you can see [examples/example.py](examples/example.py) as an example, and run with:
```sh
python3 example.py
```
