# Stripe Partitioning for the Cubed-Sphere Grid

Domain decomposition of cubed-sphere grids used in atmospheric and ocean models (CAM-SE, FV3, LFRic, E3SM, etc.).

This repository provides Fortran 90 implementations and Python wrappers for three partitioning methods:

| Method | Description |
|--------|-------------|
| **SFC** | Space-Filling Curve partitioning using composite Hilbert-Peano-Cinco curves |
| **METIS** | General-purpose multilevel graph partitioning (requires [pymetis](https://pypi.org/project/pymetis/)) |
| **Stripe** | Stripe partitioning (proposed method) — unfolds the cube panels into rectangular regions and applies optimized band decomposition |

Stripe partitioning reduces inter-process communication traffic by **15–25%** compared to SFC and by **4–6%** compared to METIS, while being deterministic and requiring no external libraries.

## Quick Start

### 1. Build Fortran code

```sh
cd f90/
make
```

Requirements: `gfortran` (or any Fortran 90 compiler).

### 2. Use in Python

**SFC partitioning**
```python
from cube_partition_sfc import CubePartitionSFC

sfc = CubePartitionSFC(ne=30, nproc=45)
nelems, cube_rank, cube_lid = sfc.make_cube_rank()
```

**Stripe partitioning**
```python
from cube_partition_stripe import CubePartitionStripe

stripe = CubePartitionStripe(ne=30, nproc=45)
nelems, cube_rank, cube_lid = stripe.make_cube_rank()
```

- `ne`: number of elements per cube edge
- `nproc`: number of processes
- `cube_rank[ei, ej, panel]`: process rank assigned to element `(ei, ej)` on `panel`

### 3. Visualize partitions

```sh
cd test/
python plot_cube_partition.py 30 45 sfc
python plot_cube_partition.py 30 45 stripe
```

### 4. Run three-way comparison (SFC vs METIS vs Stripe)

```sh
pip install pymetis   # required for METIS comparison

cd test/
python compare_three_methods.py --ne 30 --quick
python compare_three_methods.py --ne 60 --quick
python compare_three_methods.py --ne 120 --quick
```

Generate publication-quality plots:
```sh
python plot_three_methods.py
```

## Repository Structure

```
├── f90/
│   ├── cube_neighbor.f90           # Cubed-sphere neighbor connectivity
│   ├── cube_partition_sfc.f90      # SFC partitioning (Fortran)
│   ├── cube_partition_stripe.f90   # Stripe partitioning (Fortran)
│   └── makefile
├── cube_neighbor.py                # Python wrapper for cube_neighbor
├── cube_partition_sfc.py           # Python wrapper for SFC partitioning
├── cube_partition_stripe.py        # Python wrapper for stripe partitioning
├── f90wrap.py                      # Fortran-Python bridge utility
├── test/
│   ├── compare_three_methods.py    # SFC vs METIS vs Stripe comparison
│   ├── plot_three_methods.py       # Generate comparison figures
│   ├── plot_cube_partition.py      # Visualize partition shapes
│   ├── compare_cube_partition_perimeter.py
│   ├── traffic_reduction.py
│   ├── test_cube_neighbor.py
│   ├── test_cube_partition_sfc.py
│   └── test_cube_partition_stripe.py
└── README.md
```

## Method Overview

The cubed-sphere grid consists of 6 panels, each discretized into N_e x N_e elements (total: 6 N_e^2 elements). Stripe partitioning works by:

1. **Cube unfolding**: connecting the six panels into rectangular regions
2. **Band search**: a greedy algorithm determines the optimal stripe width that minimizes the mean perimeter-to-area ratio
3. **Element assignment**: elements within each stripe are assigned to processes in a zigzag pattern

An analytical lower bound for the perimeter ratio is derived and verified numerically. See the accompanying paper for details.

## Requirements

- **Fortran compiler**: gfortran (tested with GCC 12.x)
- **Python**: 3.8+
- **Python packages**: numpy, matplotlib
- **Optional**: pymetis (for METIS comparison only)

## Citation

If you use this code in your research, please cite:

```bibtex
@article{kim2026stripe,
  title={Stripe partitioning: an improved domain decomposition method for the cubed-sphere grid},
  author={Kim, Ki-Hwan and Shim, Pyoung-Seop},
  journal={Geoscientific Model Development},
  year={2026},
  note={submitted}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
