# AIMD's Block Calculation Study
Ab Initio Molecular Dynamics (AIMD) method's block calculation study.

## Installation
### Build from source
You should install [Gaussian 16](http://gaussian.com/gaussian16/) and [OpenBabel](http://openbabel.org) first. Then build AIMDBlock:
```sh
$ git clone https://github.com/njzjz/MDDatasetMaker.git
$ cd MDDatasetMaker/
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
