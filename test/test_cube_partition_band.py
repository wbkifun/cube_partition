'''

abstract : unittest of cube_partition_band.f90

history :
  2018-03-06  ki-hwan kim  start

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join
import sys
import argparse

from numpy.testing import assert_equal as equal
from numpy.testing import assert_array_equal as a_equal
from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib.patches as mp


current_dir = dirname(abspath(__file__))
lib_dir = dirname(current_dir)
sys.path.append(lib_dir)
from cube_partition_band import CubePartitionBand



def common_find_optimal_band(ne, nproc):
    start_rank = 0
    start_i = 1
    obj = CubePartitionBand(ne, nproc)
    nelems = np.ones(nproc, 'i4')*(ne*ne*6//nproc)
    if (ne*ne*6)%nproc > 0: nelems[-(ne*ne*6)%nproc:] += 1
    box = np.ones((2*ne,ne), 'i4', order='F')*(-1)
    rank, i2 = obj.find_optimal_band(start_rank, start_i, nelems, box)

    return box, rank, i2



def test_find_optimal_band_ne10():
    '''
    cube_partition_band: find_optimal_band(): ne=10, square or rectangle domain
    '''
    #
    # caution: find_optimal_band() function does not support nproc=1,2,3
    #

    ne = 10

    nproc = 4
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:ne*3//2,:], 0)
    a_equal(box[ne*3//2:,:], -1)
    equal(rank, 1)
    equal(i2, ne*3//2+1)

    nproc = 5
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:12,:], 0)
    a_equal(box[12:,:], -1)
    equal(rank, 1)
    equal(i2, 12+1)

    nproc = 6
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:ne,:], 0)
    a_equal(box[ne:,:], -1)
    equal(rank, 1)
    equal(i2, ne+1)

    nproc = 10
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:6,:], 0)
    a_equal(box[6:,:], -1)
    equal(rank, 1)
    equal(i2, 6+1)

    nproc = 12
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    '''
    a_equal(box[:ne//2,:], 0)
    a_equal(box[ne//2:,:], -1)
    equal(rank, 1)
    equal(i2, ne//2+1)
    '''
    a_equal(box[:ne,:ne//2], 0)
    a_equal(box[:ne,ne//2:], 1)
    a_equal(box[ne:,:], -1)
    equal(rank, 2)
    equal(i2, ne+1)

    nproc = 24  # 6*2*2
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//2
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:], 1)
    a_equal(box[q:,:], -1)
    equal(rank, 2)
    equal(i2, q+1)

    nproc = 150  # 6*5*5
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//5
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:q*3], 2)
    a_equal(box[:q,q*3:q*4], 3)
    a_equal(box[:q,q*4:], 4)
    equal(rank, 5)
    equal(i2, q+1)



def test_find_optimal_band_ne15():
    '''
    cube_partition_band: find_optimal_band(): ne=15, square or rectangle domain
    '''
    #
    # caution: find_optimal_band() function does not support nproc=1,2,3
    #

    ne = 15

    nproc = 6
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:ne,:], 0)
    a_equal(box[ne:,:], -1)
    equal(rank, 1)
    equal(i2, ne+1)

    nproc = 54  # 6*3*3
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//3
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:], 2)
    a_equal(box[q:,:], -1)
    equal(rank, 3)
    equal(i2, 6)

    nproc = 150  # 6*5*5
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//5
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:q*3], 2)
    a_equal(box[:q,q*3:q*4], 3)
    a_equal(box[:q,q*4:], 4)
    a_equal(box[q:,:], -1)
    equal(rank, 5)
    equal(i2, q+1)



def test_find_optimal_band_ne30():
    '''
    cube_partition_band: find_optimal_band(): ne=30, square or rectangle domain
    '''
    #
    # caution: find_optimal_band() function does not support nproc=1,2,3
    #

    ne = 30

    nproc = 4
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:ne*3//2,:], 0)
    a_equal(box[ne*3//2:,:], -1)
    equal(rank, 1)
    equal(i2, ne*3//2+1)

    nproc = 6
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:ne,:], 0)
    a_equal(box[ne:,:], -1)
    equal(rank, 1)
    equal(i2, ne+1)

    nproc = 12
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    '''
    a_equal(box[:ne//2,:], 0)
    a_equal(box[ne//2:,:], -1)
    equal(rank, 1)
    equal(i2, ne//2+1)
    '''
    a_equal(box[:ne,:ne//2], 0)
    a_equal(box[:ne,ne//2:], 1)
    a_equal(box[ne:,:], -1)
    equal(rank, 2)
    equal(i2, ne+1)

    nproc = 24  # 6*2*2
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//2
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:], 1)
    a_equal(box[q:,:], -1)
    equal(rank, 2)
    equal(i2, q+1)

    nproc = 54  # 6*3*3
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//3
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:], 2)
    a_equal(box[q:,:], -1)
    equal(rank, 3)
    equal(i2, q+1)

    nproc = 150  # 6*5*5
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//5
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:q*3], 2)
    a_equal(box[:q,q*3:q*4], 3)
    a_equal(box[:q,q*4:], 4)
    a_equal(box[q:,:], -1)
    equal(rank, 5)
    equal(i2, q+1)

    nproc = 216  # 6*6*6
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//6
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:q*3], 2)
    a_equal(box[:q,q*3:q*4], 3)
    a_equal(box[:q,q*4:q*5], 4)
    a_equal(box[:q,q*5:], 5)
    a_equal(box[q:,:], -1)
    equal(rank, 6)
    equal(i2, q+1)

    nproc = 600  # 6*10*10
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    q = ne//10
    a_equal(box[:q,:q], 0)
    a_equal(box[:q,q:q*2], 1)
    a_equal(box[:q,q*2:q*3], 2)
    a_equal(box[:q,q*3:q*4], 3)
    a_equal(box[:q,q*4:q*5], 4)
    a_equal(box[:q,q*5:q*6], 5)
    a_equal(box[:q,q*6:q*7], 6)
    a_equal(box[:q,q*7:q*8], 7)
    a_equal(box[:q,q*8:q*9], 8)
    a_equal(box[:q,q*9:], 9)
    equal(rank, 10)
    equal(i2, q+1)



def test_find_optimal_band_ne10_2():
    '''
    cube_partition_band: find_optimal_band(): ne=10, irregular domain
    '''
    ne = 10

    nproc = 7
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:8,:], 0)
    a_equal(box[8,:5], 0)
    a_equal(box[8,5:], -1)
    a_equal(box[9:,:], -1)
    equal(rank, 1)
    equal(i2, 8+1)

    nproc = 8
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:7,:], 0)
    a_equal(box[7,:5], 0)
    a_equal(box[7,5:], -1)
    a_equal(box[8:,:], -1)
    equal(rank, 1)
    equal(i2, 7+1)

    nproc = 9
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:6,:], 0)
    a_equal(box[6,:6], 0)
    a_equal(box[6,6:], -1)
    a_equal(box[7:,:], -1)
    equal(rank, 1)
    equal(i2, 6+1)

    nproc = 11
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:5,:5], 0)
    a_equal(box[5,:4], 0)
    a_equal(box[5,4:], -1)
    a_equal(box[6:,:], -1)
    equal(rank, 1)
    equal(i2, 5+1)

    nproc = 14
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:9,:4], 0)
    a_equal(box[2:8,4], 0)
    a_equal(box[:2,4], 1)
    a_equal(box[:8,5:], 1)
    a_equal(box[9,4:], -1)
    a_equal(box[10:,:], -1)
    equal(rank, 2)
    equal(i2, 9)

    nproc = 15
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:8,:5], 0)
    a_equal(box[:8,5:], 1)
    a_equal(box[8:,:], -1)
    equal(rank, 2)
    equal(i2, 8+1)

    nproc = 16
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:8,:4], 0)
    a_equal(box[2:7,4], 0)
    a_equal(box[:2,4], 1)
    a_equal(box[:7,5:], 1)
    a_equal(box[8,4:], -1)
    a_equal(box[9:,:], -1)
    equal(rank, 2)
    equal(i2, 8)

    nproc = 17
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:7,:5], 0)
    a_equal(box[:7,5:], 1)
    a_equal(box[8:,:], -1)
    equal(rank, 2)
    equal(i2, 8)

    nproc = 18
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:7,:4], 0)
    a_equal(box[2:7,4], 0)
    a_equal(box[:2,4], 1)
    a_equal(box[:7,5], 1)
    a_equal(box[6,6:], -1)
    a_equal(box[7:,:], -1)
    equal(rank, 2)
    equal(i2, 7)

    nproc = 30
    box, rank, i2 = common_find_optimal_band(ne, nproc)
    a_equal(box[:4,:5], 0)
    a_equal(box[:4,5:], 1)
    a_equal(box[4:,:], -1)
    equal(rank, 2)
    equal(i2, 4+1)



def test_make_cube_rank():
    '''
    cube_partition_band: make_cube_rank(): ne=10, square or rectangle domain
    '''
    ne = 10

    nproc = 1
    obj = CubePartitionBand(ne, nproc)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()
    a_equal(nelems, ne*ne*6)
    a_equal(cube_rank, 0)

    nproc = 2
    obj = CubePartitionBand(ne, nproc)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()
    a_equal(nelems, ne*ne*3)
    a_equal(cube_rank[:,:,-1], 0)
    a_equal(cube_rank[:,:,:2], 0)
    a_equal(cube_rank[:,:,2:5], 1)

    nproc = 3
    obj = CubePartitionBand(ne, nproc)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()
    a_equal(nelems, ne*ne*2)
    a_equal(cube_rank[:,:,-1], 0)
    a_equal(cube_rank[:,:,0], 0)
    a_equal(cube_rank[:,:,1:3], 1)
    a_equal(cube_rank[:,:,3:5], 2)

    nproc = 4
    obj = CubePartitionBand(ne, nproc)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()
    a_equal(nelems, ne*ne*6//4)
    a_equal(cube_rank[:,:,-1], 0)
    a_equal(cube_rank[:,ne//2:,0], 0)
    a_equal(cube_rank[:,:ne//2,0], 1)
    a_equal(cube_rank[:,:,1], 1)
    a_equal(cube_rank[:,:,2], 2)
    a_equal(cube_rank[:ne//2,:,3], 2)
    a_equal(cube_rank[ne//2:,:,3], 3)
    a_equal(cube_rank[:,:,4], 3)

    nproc = 6
    obj = CubePartitionBand(ne, nproc)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()
    a_equal(nelems, ne*ne)
    a_equal(cube_rank[:,:,-1], 0)
    a_equal(cube_rank[:,:,0], 1)
    a_equal(cube_rank[:,:,1], 2)
    a_equal(cube_rank[:,:,2], 3)
    a_equal(cube_rank[:,:,3], 4)
    a_equal(cube_rank[:,:,4], 5)
