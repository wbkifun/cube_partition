# Domain partitioning for cubed-sphere

This project contains Fotran codes and Python wrappers for domain partitioning a cubed-sphere grid. It supports two partitioning methods, SFC(Space Filling Curves) method and Stripe(band) method.

## Quick start

### 1. Build Fortran codes
```sh
$ cd f90/
$ make
```

### 2. Get partitioning results in Python
**Using SFC method**
```python
from cube_partition_sfc import CubePartitionSFC

sfc = CubePartitionSFC(ne, nproc)
nelems, cube_rank, cube_lid = sfc.make_cube_rank()
```

**Using Stripe(band) method**
```python
from cube_partition_band import CubePartitionBand

band = CubePartitionBand(ne, nproc)
nelems, cube_rank, cube_lid = band.make_cube_rank()
```

### 3. Plot partitioned cubed-sphere grid
```sh
$ cd test/
$ python plot_cube_partition.py -h
usage: plot_cube_partition.py [-h] [--save] [--rank_fontsize RANK_FONTSIZE] ne nproc {sfc,band}

positional arguments:
  ne                    number of elements
  nproc                 number of processors
  {sfc,band}            partitioning method

optional arguments:
  -h, --help            show this help message and exit
  --save                save as png format (default: False)
  --rank_fontsize RANK_FONTSIZE
                        fontsize of rank numbers (default: 0)

$ python plot_cube_partition.py 10 30 sfc
$ python plot_cube_partition.py 10 30 band
```
