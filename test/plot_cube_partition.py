'''

abstract : plot the partitioning of the cubed-sphere

history :
  2018-03-29  ki-hwan kim  split from test_cube_partition_band.py

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join
import sys
import argparse

from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib.patches as mp


current_dir = dirname(abspath(__file__))
lib_dir = dirname(current_dir)
sys.path.append(lib_dir)
from cube_partition_sfc import CubePartitionSFC
from cube_partition_band import CubePartitionBand



def discrete_cmap(n, base_cmap=None):
    ''' 
    create an N-bin discrete colormap from the specified input map
    insert white at the front
    '''

    base = plt.cm.get_cmap(base_cmap)
    color_list = np.zeros((n,4), 'f4')
    color_list[0,:] = [1,1,1,1]  # white
    color_list[1:,:] = base(np.linspace(0,1,n-1))
    cmap_name = base.name + str(n)
    return base.from_list(cmap_name, color_list, n)



def discrete7_cmap(n):
    ''' 
    create an 7-bin discrete colormap
    '''

    color_list = np.zeros((7,3), 'f4')
    color_list[0,:] = [1,1,1]  # white
    color_list[1,:] = [1,0,0]  # red
    color_list[2,:] = [0,1,0]  # green
    color_list[3,:] = [0,0,1]  # blue
    color_list[4,:] = [1,1,0]  # yellow
    color_list[5,:] = [1,0,1]  # magenta
    color_list[6,:] = [0,1,1]  # cyan

    if n <= 6:
        return ListedColormap(color_list[:n+1,:])
    else:
        color_list2 = np.zeros((n+1,3), 'f4')
        color_list2[:7,:] = color_list
        color_list2[7:,:] = np.random.rand(n-6,3)
        return ListedColormap(color_list2)



def plot_find_optimal_band():
    '''
    cube_partition_band: find_optimal_band(): plot the box
    '''
    ne, nproc = 10, 37
    obj = CubePartitionBand(ne, nproc)

    box = np.ones((2*ne,ne), 'i4', order='F')*(-1)
    start_rank = 0
    start_i = 1
    nelems = np.ones(nproc, 'i4')*(ne*ne*6//nproc)
    if (ne*ne*6)%nproc > 0: nelems[-(ne*ne*6)%nproc:] += 1
    rank, i2 = obj.find_optimal_band(start_rank, start_i, nelems, box)
    print('nelems', nelems)
    print('rank={}, i2={}'.format(rank, i2))

    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1,1,1)
    mycmap = discrete_cmap(nproc+1, 'prism')  # Paired, prism, Set1
    ax.imshow(box.T, origin='lower', cmap=mycmap, interpolation='nearest')
    ax.set_xticks(np.arange(0.5,2*ne-0.5))
    ax.set_yticks(np.arange(0.5,ne-0.5))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(ls='-')
    plt.tight_layout()
    plt.show()



def check_perimeter_ratio():
    '''
    cube_partition_band: check_perimeter_ratio()
    '''
    opt = 1

    if opt == 1:
        #
        # plot the perimeter ratios
        #
        ne = 191
        nproc = ne  # maximum number of processes
        nelem = 53

        obj = CubePartitionBand(ne, nproc)
        nelems = np.ones(ne, 'i4')*nelem
        ratios = np.zeros(ne, 'f8')
        box = np.ones((2*ne, ne), 'i4', order='F')*(-1)

        for end_rank in range(ne):
            i12 = np.array([1, 1, ne, 0], 'i4')  # (i1, i2, band_elem, remain_elem)
            box[:,:] = -1
            ratio = obj.calc_perimeter_ratio(0, end_rank, nelems, i12, box)
            ratios[end_rank] = ratio
            print('block={}, {}'.format(end_rank+1, ratio))

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(1,1,1)
        ax.plot(np.arange(1,ne+1), ratios, '-o')
        ax.set_ylim(0.4,1)
        plt.tight_layout()
        plt.show()

    elif opt == 2:
        #
        # plot the box
        #
        ne = 10
        nproc = ne  # maximum number of processes
        nelem = 13

        obj = CubePartitionBand(ne, nproc)
        nelems = np.ones(ne, 'i4')*nelem
        i12 = np.array([1, 1, ne, 0], 'i4')  # (i1, i2, band_elem, remain_elem)
        box = np.ones((2*ne, ne), 'i4', order='F')*(-1)

        end_rank = 2
        ratio = obj.calc_perimeter_ratio(0, end_rank, nelems, i12, box)
        print('ne={}, block={}, ratio={}'.format(ne, end_rank+1, ratio))

        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(1,1,1)
        mycmap = discrete_cmap(nproc+1, 'prism') # Paired, prism, Set1
        ax.imshow(box.T, origin='lower', cmap=mycmap, interpolation='nearest')

        ax.set_yticks(np.arange(0.5,ne-0.5))
        ax.set_xticks(np.arange(0.5,2*ne-0.5))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(ls='-')
        plt.tight_layout()
        plt.show()



def fit_perimeter_ratio():
    '''
    fitting the perimeter ratios
    '''

    #
    # analytic formula
    #
    w = 200  # width
    s = 16  # block size

    r = np.zeros(w, 'f8')
    for x in range(1,w+1):
        xp = x/w  # normalized x'=x/w
        if s*xp <= 1: xp = 1/s
        r[x-1] = 2*(1/(s*xp) + xp)

    #
    # experiments
    #
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)

    for ne in [240, 120, 60]:
        nproc = ne  # maximum number of processes
        nelem = s

        obj = CubePartitionBand(ne, nproc)
        nelems = np.ones(ne, 'i4')*nelem
        ratios = np.zeros(ne, 'f8')
        box = np.ones((2*ne, ne), 'i4', order='F')*(-1)

        for end_rank in range(ne):
            i12 = np.array([1, 1, ne, 0], 'i4')  # (i1, i2, band_elem, remain_elem)
            box[:,:] = -1
            ratio = obj.calc_perimeter_ratio(0, end_rank, nelems, i12, box)
            ratios[end_rank] = ratio
        #print('ne, ratios', ne, ratios)

        xx = np.arange(1,ne+1)/ne
        ax.plot(xx, ratios, 'o', label=r'$N_e={}$'.format(ne))

    xx = np.arange(1,w+1)/w
    ax.plot(xx, r, '-k', lw=2, label='Lower bound')
    #ax.plot([max_xp,max_xp], [0,1], ':k', lw=2)

    ax.set_ylim(0.9, 2.2)
    ax.set_xlabel(r'$N_{block}/N_{e}$')
    ax.set_ylabel('Perimeter/Area')
    ax.legend(loc='lower right', numpoints=1)
    plt.tight_layout()
    plt.savefig('png/minimum_perimeter_ratio.png', dpi=300)
    plt.show()




def plot_cube_partition(ne, nproc, method, save, rank_fontsize):
    '''
    cube_partition_band: plot the cube_rank array
    '''
    obj = CubePartitionBand(ne, nproc)

    if method == 'sfc':
        sfc = CubePartitionSFC(ne, nproc)
        nelems, cube_rank, cube_lid = sfc.make_cube_rank()
    elif method == 'band':
        nelems, cube_rank, cube_lid = obj.make_cube_rank()

    cube_color = obj.make_cube_color(cube_rank)

    perimeter_ratio, num_nbrs = obj.global_perimeter_ratio(cube_rank)
    print('ne={}, nproc={}, perimter_ratio={}'.format(ne, nproc, perimeter_ratio))

    box = np.zeros((4*ne,3*ne), 'i4', order='F')
    box[:ne,2*ne:] = cube_color[:,:,5]
    box[:ne,ne:2*ne] = cube_color[:,:,0]
    box[ne:2*ne,ne:2*ne] = cube_color[:,:,1]
    box[2*ne:3*ne,ne:2*ne] = cube_color[:,:,2]
    box[3*ne:,ne:2*ne] = cube_color[:,:,3]
    box[3*ne:,:ne] = np.rot90(cube_color[:,:,4], 3)

    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1,1,1)
    #mycmap = discrete_cmap(nproc+1, 'prism') # Paired, prism, Set1
    mycmap = discrete7_cmap(box.max())
    ax.imshow(box.T, origin='lower', cmap=mycmap, interpolation='nearest')
    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.axis('off')

    #
    # draw element lines
    #
    vlines = np.zeros((4*ne-1,2,2), 'f4')
    hlines = np.zeros((3*ne-1,2,2), 'f4')

    for i in range(4*ne-1):
        x = i+0.5
        y1, y2 = -0.5, 3*ne+0.5

        if x < 3*ne: y1 = ne-0.5
        if x > ne-1: y2 = 2*ne-0.5

        vlines[i,:,:] = [(x, y1), (x, y2)]

    for j in range(3*ne-1):
        x1, x2 = -0.5, 4*ne+0.5
        y = j+0.5

        if y < ne: x1 = 3*ne-0.5
        if y > 2*ne-1: x2 = ne-0.5

        hlines[j,:,:] = [(x1, y), (x2, y)]

    lcv = mc.LineCollection(vlines, colors='k')
    lch = mc.LineCollection(hlines, colors='k')
    ax.add_collection(lcv)
    ax.add_collection(lch)

    #
    # draw panel boxes
    #
    kwds = {'fill':False, 'linewidth':2}
    rects = [mp.Rectangle((-0.5, 2*ne-0.5), ne, ne, **kwds),
             mp.Rectangle((-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((ne-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((2*ne-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((3*ne-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((3*ne-0.5, -0.5), ne, ne, **kwds)]
    for rect in rects: ax.add_patch(rect)

    #
    # draw rank numbers
    #
    if rank_fontsize != 0:
        box2 = np.ones((4*ne,3*ne), 'i4', order='F')*(-1)
        box2[:ne,2*ne:] = cube_rank[:,:,5]
        box2[:ne,ne:2*ne] = cube_rank[:,:,0]
        box2[ne:2*ne,ne:2*ne] = cube_rank[:,:,1]
        box2[2*ne:3*ne,ne:2*ne] = cube_rank[:,:,2]
        box2[3*ne:,ne:2*ne] = cube_rank[:,:,3]
        box2[3*ne:,:ne] = np.rot90(cube_rank[:,:,4], 3)
        for (j, i), label in np.ndenumerate(box2.T):
            if label >= 0:
                ax.text(i, j, label, ha='center', va='center', fontsize=rank_fontsize)

    #
    # plot and save
    #
    plt.tight_layout()
    if save:
        if rank_fontsize !=0:
            fname = 'cube_partition.ne{}_nproc{}.{}.number.png'.format(ne, nproc, method)
        else:
            fname = 'cube_partition.ne{}_nproc{}.{}.png'.format(ne, nproc, method)
        plt.savefig('png/'+fname, dpi=300)

    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--save', action='store_true', help='save as png format')
    parser.add_argument('--rank_fontsize', type=int, default=0, help='fontsize of rank numbers')
    parser.add_argument('ne', type=int, help='number of elements')
    parser.add_argument('nproc', type=int, help='number of processors')
    parser.add_argument('method', type=str, choices=['sfc','band'], help='partitioning method')
    args = parser.parse_args()

    plot_cube_partition(args.ne, args.nproc, args.method, args.save, args.rank_fontsize)
