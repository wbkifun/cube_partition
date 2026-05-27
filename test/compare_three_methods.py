'''
Compare three partitioning methods on the cubed-sphere:
  1. Space-Filling Curve (SFC)
  2. Stripe (Band) partitioning
  3. METIS graph partitioning

Metrics:
  - Mean perimeter ratio (boundary elements / total elements per partition)
  - Communication ratio for Spectral Element Method (Np=4)
  - Total communication traffic
  - Partitioning wall-clock time
'''

from os.path import dirname, abspath, join
import sys
import time
import numpy as np
import pymetis

current_dir = dirname(abspath(__file__))
lib_dir = dirname(current_dir)
sys.path.append(lib_dir)

from cube_partition_sfc import CubePartitionSFC
from cube_partition_stripe import CubePartitionStripe
from cube_neighbor import CubeNeighbor


def build_cubed_sphere_adjacency(ne):
    '''
    Build the adjacency list for the cubed-sphere grid.
    Each element (ei, ej, panel) is a node.
    Edges connect 4-neighbors (N, S, E, W) across panel boundaries.
    Returns adjacency list suitable for pymetis.
    '''
    cn = CubeNeighbor()
    n_total = 6 * ne * ne

    def eid(ei, ej, p):
        return (p - 1) * ne * ne + (ej - 1) * ne + (ei - 1)

    adjacency = [[] for _ in range(n_total)]

    for p in range(1, 7):
        for ej in range(1, ne + 1):
            for ei in range(1, ne + 1):
                me = eid(ei, ej, p)
                for (di, dj) in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nei, nej, np_ei, np_ej = ei + di, ej + dj, 0, 0
                    ret = cn.convert_nbr_eij(ne, nei, nej, p)
                    nbr_ei, nbr_ej, nbr_p = int(ret[0]), int(ret[1]), int(ret[2])
                    if nbr_p >= 1:
                        nbr = eid(nbr_ei, nbr_ej, nbr_p)
                        if nbr not in adjacency[me]:
                            adjacency[me].append(nbr)

    return adjacency


def metis_partition(ne, nproc, adjacency):
    '''
    Partition the cubed-sphere using METIS.
    Returns cube_rank array (ne, ne, 6) in Fortran order.
    '''
    n_total = 6 * ne * ne

    if nproc == 1:
        membership = [0] * n_total
    else:
        _, membership = pymetis.part_graph(nproc, adjacency=adjacency)

    cube_rank = np.zeros((ne, ne, 6), 'i4', order='F')
    for p in range(6):
        for ej in range(ne):
            for ei in range(ne):
                idx = p * ne * ne + ej * ne + ei
                cube_rank[ei, ej, p] = membership[idx]

    return cube_rank


def compute_metrics(ne, nproc, cube_rank, band_obj, ngq=4):
    '''
    Compute perimeter ratio and communication ratio for a given partition.
    '''
    pr, num_nbrs = band_obj.global_perimeter_ratio(cube_rank)
    cr, num_pts = band_obj.global_communication_ratio(ngq, cube_rank)

    ratios = num_nbrs[1, :] / np.maximum(num_nbrs[0, :], 1)
    pr_std = ratios.std()

    comm_ratios = num_pts[1, :] / np.maximum(num_pts[0, :], 1)
    cr_std = comm_ratios.std()

    total_comm = num_pts[1, :].sum()

    return {
        'perimeter_ratio_mean': pr,
        'perimeter_ratio_std': pr_std,
        'comm_ratio_mean': cr,
        'comm_ratio_std': cr_std,
        'total_comm': total_comm,
    }


def run_comparison(ne, nproc_list, ngq=4):
    '''
    Run the three-way comparison for a given Ne and list of Nproc values.
    '''
    print(f'\n{"="*70}')
    print(f'  Cubed-Sphere Partitioning Comparison: Ne={ne}, Ngq={ngq}')
    print(f'  Total elements: {6*ne*ne}')
    print(f'{"="*70}')

    print('\nBuilding adjacency list for METIS...')
    t0 = time.time()
    adjacency = build_cubed_sphere_adjacency(ne)
    t_adj = time.time() - t0
    print(f'  Adjacency built in {t_adj:.2f}s')

    header = (f'{"Nproc":>6} | {"Method":>8} | {"P_mean":>8} | {"P_std":>8} | '
              f'{"CR_mean":>8} | {"CR_std":>8} | {"TotalComm":>10} | {"Time(ms)":>8}')
    print(f'\n{header}')
    print('-' * len(header))

    results = []

    for nproc in nproc_list:
        if nproc < 1:
            continue

        band_obj = CubePartitionStripe(ne, nproc)

        # --- SFC ---
        t0 = time.time()
        sfc_obj = CubePartitionSFC(ne, nproc)
        _, cr_sfc, _ = sfc_obj.make_cube_rank()
        t_sfc = (time.time() - t0) * 1000

        m_sfc = compute_metrics(ne, nproc, cr_sfc, band_obj, ngq)
        m_sfc['time_ms'] = t_sfc
        m_sfc['method'] = 'SFC'
        m_sfc['nproc'] = nproc
        results.append(m_sfc)

        print(f'{nproc:>6} | {"SFC":>8} | {m_sfc["perimeter_ratio_mean"]:>8.4f} | '
              f'{m_sfc["perimeter_ratio_std"]:>8.4f} | {m_sfc["comm_ratio_mean"]:>8.4f} | '
              f'{m_sfc["comm_ratio_std"]:>8.4f} | {m_sfc["total_comm"]:>10} | {t_sfc:>8.1f}')

        # --- Stripe ---
        t0 = time.time()
        _, cr_band, _ = band_obj.make_cube_rank()
        t_band = (time.time() - t0) * 1000

        m_band = compute_metrics(ne, nproc, cr_band, band_obj, ngq)
        m_band['time_ms'] = t_band
        m_band['method'] = 'Stripe'
        m_band['nproc'] = nproc
        results.append(m_band)

        print(f'{"":>6} | {"Stripe":>8} | {m_band["perimeter_ratio_mean"]:>8.4f} | '
              f'{m_band["perimeter_ratio_std"]:>8.4f} | {m_band["comm_ratio_mean"]:>8.4f} | '
              f'{m_band["comm_ratio_std"]:>8.4f} | {m_band["total_comm"]:>10} | {t_band:>8.1f}')

        # --- METIS ---
        t0 = time.time()
        cr_metis = metis_partition(ne, nproc, adjacency)
        t_metis = (time.time() - t0) * 1000

        m_metis = compute_metrics(ne, nproc, cr_metis, band_obj, ngq)
        m_metis['time_ms'] = t_metis
        m_metis['method'] = 'METIS'
        m_metis['nproc'] = nproc
        results.append(m_metis)

        print(f'{"":>6} | {"METIS":>8} | {m_metis["perimeter_ratio_mean"]:>8.4f} | '
              f'{m_metis["perimeter_ratio_std"]:>8.4f} | {m_metis["comm_ratio_mean"]:>8.4f} | '
              f'{m_metis["comm_ratio_std"]:>8.4f} | {m_metis["total_comm"]:>10} | {t_metis:>8.1f}')

        # --- Summary for this nproc ---
        sfc_tc = m_sfc['total_comm']
        band_tc = m_band['total_comm']
        metis_tc = m_metis['total_comm']
        print(f'{"":>6} | {"":>8} | Traffic reduction vs SFC:  '
              f'Stripe={((sfc_tc-band_tc)/sfc_tc*100):>5.1f}%  '
              f'METIS={((sfc_tc-metis_tc)/sfc_tc*100):>5.1f}%')
        print('-' * len(header))

    return results


def save_results(results, ne, fname=None):
    '''Save results to numpy file.'''
    if fname is None:
        fname = f'data/comparison_three_methods_ne{ne}.npy'
    np.save(fname, results)
    print(f'\nResults saved to {fname}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ne', type=int, default=30, help='number of elements per edge')
    parser.add_argument('--ngq', type=int, default=4, help='number of Gauss quadrature points')
    parser.add_argument('--quick', action='store_true', help='quick test with fewer nproc values')
    args = parser.parse_args()

    ne = args.ne

    if args.quick:
        nproc_list = [6, 12, 24, 36, 48, 54, 72, 96]
    else:
        total_elems = 6 * ne * ne
        min_elems_per_proc = 4
        max_nproc = min(total_elems // min_elems_per_proc, 200)
        if ne <= 30:
            step = max(1, max_nproc // 20)
            nproc_list = list(range(step, max_nproc + 1, step))
        else:
            step = max(1, max_nproc // 15)
            nproc_list = list(range(step, max_nproc + 1, step))

    nproc_list = [n for n in nproc_list if n >= 4]

    results = run_comparison(ne, nproc_list, args.ngq)
    save_results(results, ne)
