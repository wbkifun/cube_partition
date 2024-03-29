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



def get_communication_ratios(ne, ngq, max_nproc):
    '''
    save the communication/computation ratios along the number of processes
    to a numpy binary file
    '''
    fpath_base = 'data/comm_ratios_ne{}np{}'.format(ne, ngq)
    fpath_sfc = 'data/comm_ratios_ne{}np{}_sfc.npy'.format(ne, ngq)
    fpath_band = 'data/comm_ratios_ne{}np{}_band.npy'.format(ne, ngq)
    fpath_total = 'data/comm_total_ne{}np{}.npy'.format(ne, ngq)

    if exists(fpath_sfc) and exists(fpath_band) and exists(fpath_total):
        return np.load(fpath_sfc), np.load(fpath_band), np.load(fpath_total)

    sfc_std = np.zeros((max_nproc,3), 'f4')  # (mean, mean-std, mean+std)
    band_std = np.zeros((max_nproc,3), 'f4')
    total_comm = np.zeros((max_nproc,2), 'f4')

    for i in range(max_nproc):
        nproc = i + 1

        sfc = CubePartitionSFC(ne, nproc)
        band = CubePartitionBand(ne, nproc)

        nelems1, cube_rank1, cube_lid1 = sfc.make_cube_rank()
        nelems2, cube_rank2, cube_lid2 = band.make_cube_rank()

        mean_sfc, num_pts_sfc = band.global_communication_ratio(ngq, cube_rank1)
        mean_band, num_pts_band = band.global_communication_ratio(ngq, cube_rank2)
        print('{}\t{}\t{}'.format(nproc, mean_sfc, mean_band))

        ratios_sfc = num_pts_sfc[1,:]/num_pts_sfc[0,:]
        ratios_band = num_pts_band[1,:]/num_pts_band[0,:]

        # standrad deviation
        std1 = ratios_sfc.std()
        std2 = ratios_band.std()
        sfc_std[i,:] = (mean_sfc, mean_sfc-std1, mean_sfc+std1)
        band_std[i,:] = (mean_band, mean_band-std2, mean_band+std2)

        # total communication count
        total_comm[i,0] = sum(num_pts_sfc[1,:])
        total_comm[i,1] = sum(num_pts_band[1,:])

    np.save(fpath_sfc, sfc_std)
    np.save(fpath_band, band_std)
    np.save(fpath_total, total_comm)

    return sfc_std, band_std, total_comm



def compare_communication_ratios(ne, ngq, save):
    '''
    compare the communication/computation ratios between two partitioing methods
    '''
    #xmax = 1/16  # max nproc ratio, 1/16
    xmax = 1/32  # max nproc ratio, 1/32

    max_nproc = int(6*ne*ne*xmax)
    gap_nproc = 6*(ne//60)**2
    nx = max_nproc//gap_nproc
    print('ne={}, np={}, max_nproc={}, nx={}'.format(ne, ngq, max_nproc, nx))

    sfc_ratios, band_ratios, total_comm = get_communication_ratios(ne, ngq, max_nproc)

    sfc_plot = np.zeros((nx,3), 'f4')  # (mean, min, max)
    band_plot = np.zeros((nx,3), 'f4')  # (mean, min, max)

    for i in range(nx):
        nproc = (i+1)*gap_nproc
        sfc_plot[i,:] = sfc_ratios[nproc-1,:]
        band_plot[i,:] = band_ratios[nproc-1,:]


    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim(gap_nproc, max_nproc)

    x = np.arange(gap_nproc,max_nproc+1,gap_nproc)

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
    ax.set_ylabel('Communication/Computation')
    ax.legend([l3, l2], 
              ['Partitioning with Space-Filling Curves',
               'Stripe partitioning'],
              loc='lower right',
              fontsize='medium')

    plt.tight_layout()
    if save:
        fname = 'compare_sfc_band_ne{}np{}_maxnproc{}_comm.png'.format(ne, ngq, max_nproc)
        plt.savefig('png/'+fname, dpi=300)
    plt.show()



def compare_communication_ratios2(ne_list, ngq, save):
    '''
    compare the communication/computation ratios between two partitioing methods
    '''

    max_nproc = 500

    total_plots = np.zeros((len(ne_list), max_nproc), 'f4')

    for i, ne in enumerate(ne_list):
        sfc_ratios, band_ratios, total_comm = get_communication_ratios(ne, ngq, max_nproc)
        for j in range(max_nproc):
            total_plots[i,j] = (total_comm[j,0] - total_comm[j,1])/total_comm[j,0]*100


    fig = plt.figure(figsize=(7,4))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim(2, max_nproc)

    for i, ne in enumerate(ne_list):
        x = np.arange(2, max_nproc+1)
        ax.plot(x, total_plots[i,1:], 'o-', lw=1, ms=2, label=r'$N_{e}=%d$'%(ne))

    #ax.set_ylim(0, 1)
    ax.set_xlabel('Number of processes')
    ax.set_ylabel('Total traffic reduction [%]')
    ax.set_title(r'Spectral Element Method ($N_p=4$)')
    ax.legend(loc='lower center', fontsize='medium')

    plt.tight_layout()
    if save:
        fname = 'traffic_reduction_np{}_maxnproc{}.png'.format(ngq, max_nproc)
        plt.savefig('png/'+fname, dpi=300)
    plt.show()



if __name__ == '__main__':
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ne', type=int, help='number of elements')
    parser.add_argument('np', type=int, nargs='?', default=4, help='number of quadrature points')
    parser.add_argument('--save', action='store_true', help='save as png format')
    args = parser.parse_args()

    compare_communication_ratios(args.ne, args.np, args.save)
    '''
    compare_communication_ratios2([60, 120, 240], 4, True)
