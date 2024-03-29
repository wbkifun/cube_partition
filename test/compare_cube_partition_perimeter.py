'''

abstract : compare two partitioning methods

history :
  2018-03-23  ki-hwan kim  start

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join, exists
import sys
import argparse

from numpy.testing import assert_equal as equal
from numpy.testing import assert_array_equal as a_equal
import numpy as np
import matplotlib.pyplot as plt


current_dir = dirname(abspath(__file__))
lib_dir = dirname(current_dir)
sys.path.append(lib_dir)
from cube_partition_band import CubePartitionBand
from cube_partition_sfc import CubePartitionSFC



def get_perimeter_ratios(ne, max_nproc, spread):
    '''
    save the perimeter ratios along the number of processes
    to a numpy binary file
    '''
    fpath_base = 'data/perimeter_ratios_ne{}'.format(ne)
    fpath_sfc = fpath_base + '_sfc_{}.npy'.format(spread)
    fpath_band = fpath_base + '_band_{}.npy'.format(spread)

    if exists(fpath_sfc) and exists(fpath_band):
        return np.load(fpath_sfc), np.load(fpath_band)

    sfc_minmax = np.zeros((max_nproc,3), 'f4') # (mean, min, max)
    sfc_std = np.zeros((max_nproc,3), 'f4')    # (mean, mean-std, mean+std)
    band_minmax = np.zeros((max_nproc,3), 'f4')
    band_std = np.zeros((max_nproc,3), 'f4')

    for i in range(max_nproc):
        nproc = i + 1

        sfc = CubePartitionSFC(ne, nproc)
        band = CubePartitionBand(ne, nproc)

        nelems1, cube_rank1, cube_lid1 = sfc.make_cube_rank()
        nelems2, cube_rank2, cube_lid2 = band.make_cube_rank()

        mean_sfc, num_nbrs_sfc = band.global_perimeter_ratio(cube_rank1)
        mean_band, num_nbrs_band = band.global_perimeter_ratio(cube_rank2)
        print('{}\t{}\t{}'.format(nproc, mean_sfc, mean_band))

        ratios_sfc = num_nbrs_sfc[1,:]/num_nbrs_sfc[0,:]
        ratios_band = num_nbrs_band[1,:]/num_nbrs_band[0,:]

        # minmax
        sfc_minmax[i,:] = (mean_sfc, ratios_sfc.min(), ratios_sfc.max())
        band_minmax[i,:] = (mean_band, ratios_band.min(), ratios_band.max())

        # standrad deviation
        std1 = ratios_sfc.std()
        std2 = ratios_band.std()
        sfc_std[i,:] = (mean_sfc, mean_sfc-std1, mean_sfc+std1)
        band_std[i,:] = (mean_band, mean_band-std2, mean_band+std2)

    fpath = 'data/perimeter_ratios_ne{}_sfc_minmax.npy'.format(ne)
    np.save(fpath, sfc_minmax)
    fpath = 'data/perimeter_ratios_ne{}_sfc_std.npy'.format(ne)
    np.save(fpath, sfc_std)
    fpath = 'data/perimeter_ratios_ne{}_band_minmax.npy'.format(ne)
    np.save(fpath, band_minmax)
    fpath = 'data/perimeter_ratios_ne{}_band_std.npy'.format(ne)
    np.save(fpath, band_std)

    if spread == 'minmax':
        return sfc_minmax, band_minmax
    elif spread == 'std':
        return sfc_std, band_std



def compare_perimeter_ratios(ne, spread, save):
    '''
    compare the perimeter ratios between two partitioing methods
    spread: 'std', 'minmax'
    '''
    assert spread in ['std', 'minmax']

    xmax = 1/16  # max nproc ratio, 1/16
    #xmax = 1/32  # max nproc ratio, 1/32

    max_nproc = int(6*ne*ne*xmax)
    gap_nproc = 6*(ne//60)**2
    nx = max_nproc//gap_nproc
    print('ne={}, max_nproc={}, nx={}'.format(ne, max_nproc, nx))

    sfc_ratios, band_ratios = get_perimeter_ratios(ne, max_nproc, spread)

    low_limit = np.zeros(nx, 'f4')
    #low_limit2 = np.zeros(nx, 'f4')
    sfc_plot = np.zeros((nx,3), 'f4')  # (mean, min, max)
    band_plot = np.zeros((nx,3), 'f4')  # (mean, min, max)

    for i in range(nx):
        nproc = (i+1)*gap_nproc

        sfc_plot[i,:] = sfc_ratios[nproc-1,:]
        band_plot[i,:] = band_ratios[nproc-1,:]

        area = (ne*ne*6)/nproc
        low_limit[i] = 4*np.sqrt(area)/area
        #low_limit2[i] = 2*np.sqrt(3*area)/area


    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim(gap_nproc, max_nproc)

    x = np.arange(gap_nproc,max_nproc+1,gap_nproc)

    l1, = ax.plot(x, low_limit, 'k-', lw=2)
    #ax.plot(x, low_limit2, 'k:' )

    l2, = ax.plot(x, band_plot[:,0], 'ro-', ms=3)
    ax.fill_between(x, band_plot[:,1], band_plot[:,2], facecolor='r', alpha=0.3)

    l3, = ax.plot(x, sfc_plot[:,0], 'bo-', ms=3)
    ax.fill_between(x, sfc_plot[:,1], sfc_plot[:,2], facecolor='b', alpha=0.3)

    if ne in [60, 120, 240]:
        seq = {60:'a', 120:'b', 240:'c'}
        ax.set_title(r'({}) $Ne=${}'.format(seq[ne], ne))
    else:
        ax.set_title('ne={}'.format(ne))

    #ax.set_ylim(0, 1)
    ax.set_xlabel('Number of processes')
    ax.set_ylabel('Perimeter/Area')
    ax.legend([l3, l2, l1], 
              ['Partitioning with Space-Filling Curves',
               'Stripe partitioning',
               'Lower limit of stripe partitioning'],
              loc='lower right',
              fontsize='medium')

    plt.tight_layout()
    if save:
        fname = 'compare_sfc_band_ne{}_maxnproc{}.png'.format(ne, max_nproc)
        plt.savefig('png/'+fname, dpi=300)
    plt.show()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ne', type=int, help='number of elements')
    parser.add_argument('--spread', type=str, default='std', choices=['std', 'minmax'], help='speread type')
    parser.add_argument('--save', action='store_true', help='save as png format')
    args = parser.parse_args()

    compare_perimeter_ratios(args.ne, args.spread, args.save)
