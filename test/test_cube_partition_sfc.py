'''

abstract : unittest of cube_partition_sfc.f90

history :
  2017-04-06  ki-hwan kim  initial setup
  2017-06-16  ki-hwan kim  add make_elem_coord()

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
from cube_partition_sfc import CubePartitionSFC
from cube_partition_band import CubePartitionBand



def test_rot():
    '''
    cube_partition_sfc: rot()
    '''
    obj = CubePartitionSFC(ne=30, nproc=10)

    arr = np.array([[1,2],[3,4]], dtype='i4', order='F')
    a_equal(obj.rot(1, arr), np.rot90(arr))
    a_equal(obj.rot(2, arr), np.rot90(np.rot90(arr)))
    a_equal(obj.rot(3, arr), np.rot90(np.rot90(np.rot90(arr))))

    arr = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype='i4', order='F')
    a_equal(obj.rot(1, arr), np.rot90(arr))
    a_equal(obj.rot(2, arr), np.rot90(np.rot90(arr)))
    a_equal(obj.rot(3, arr), np.rot90(np.rot90(np.rot90(arr))))



def test_inv_x():
    '''
    cube_partition_sfc: inv_x() and inv_y()
    '''
    obj = CubePartitionSFC(ne=30, nproc=10)

    arr = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype='i4', order='F')
    a_equal(obj.inv_x(arr), arr[::-1,:])
    a_equal(obj.inv_y(arr), arr[:,::-1])



def test_make_sfcs():
    '''
    cube_partition_sfc: make_sfcs()
    '''
    obj = CubePartitionSFC(ne=30, nproc=10)
    hilbert, peano, cinco = obj.make_sfcs()
    a_equal(hilbert[:,:,0], np.array([[1,2],[4,3]]))
    a_equal(hilbert[:,:,1], np.array([[1,4],[2,3]]))
    a_equal(hilbert[:,:,2], np.array([[3,2],[4,1]]))
    a_equal(hilbert[:,:,3], np.array([[3,4],[2,1]]))



def test_find_factors():
    '''
    cube_partition_sfc: find_factors()
    '''
    obj = CubePartitionSFC(ne=2*3*5, nproc=1)
    a_equal(obj.find_factors(), [2,3,5])

    '''
    obj = CubePartitionSFC(ne=2*3*5, nproc=6*2*2)
    a_equal(obj.find_factors(), [3,5,2])

    obj = CubePartitionSFC(ne=2*3*5, nproc=6*3*3)
    a_equal(obj.find_factors(), [2,5,3])
    '''

    obj = CubePartitionSFC(ne=120, nproc=10)
    a_equal(obj.find_factors(), [2,2,2,3,5])

    '''
    obj = CubePartitionSFC(ne=120, nproc=6*2*2*2)
    a_equal(obj.find_factors(), [3,5,2,2,2])
    '''



def test_make_panel_sfc():
    '''
    cube_partition_sfc: make_panel_sfc()
    '''
    obj = CubePartitionSFC(ne=2*2, nproc=1)
    expect = np.array([ \
            [ 1, 4, 5, 6],
            [ 2, 3, 8, 7],
            [15,14, 9,10], 
            [16,13,12,11]])
    a_equal(obj.make_panel_sfc(), expect)

    obj = CubePartitionSFC(ne=2*3, nproc=1)
    expect = np.array([ \
            [ 1, 4, 5, 8, 9,10],
            [ 2, 3, 6, 7,12,11],
            [31,30,27,26,13,14],
            [32,29,28,25,16,15],
            [33,34,23,24,17,18],
            [36,35,22,21,20,19]]) 
    a_equal(obj.make_panel_sfc(), expect)

    obj = CubePartitionSFC(ne=2*5, nproc=1)
    expect = np.array([ \
            [ 1, 2,31,32,33,36,37,40,41,42],
            [ 4, 3,30,29,34,35,38,39,44,43],
            [ 5, 6,27,28,23,22,51,50,45,46],
            [ 8, 7,26,25,24,21,52,49,48,47],
            [ 9,12,13,16,17,20,53,56,57,58],
            [10,11,14,15,18,19,54,55,60,59],
            [95,94,91,90,79,78,75,74,61,62],
            [96,93,92,89,80,77,76,73,64,63],
            [97,98,87,88,81,82,71,72,65,66],
           [100,99,86,85,84,83,70,69,68,67]]) 
    a_equal(obj.make_panel_sfc(), expect)



def test_make_global_sfc_2():
    '''
    cube_partition_sfc: make_global_sfc(ne=2)
    '''
    expect_p1 = np.array( \
            [[ 2, 1], \
             [ 3, 4]])

    expect_p2 = np.array( \
            [[ 6, 5], \
             [ 7, 8]])

    expect_p3 = np.array( \
            [[21,22], \
             [24,23]])

    expect_p4 = np.array( \
            [[16,13], \
             [15,14]])

    expect_p5 = np.array( \
            [[17,18], \
             [20,19]])

    expect_p6 = np.array( \
            [[11,12], \
             [10, 9]])

    obj = CubePartitionSFC(ne=2, nproc=1)
    cube_gid = obj.make_global_sfc()
    a_equal(cube_gid[:,:,0], expect_p1)
    a_equal(cube_gid[:,:,1], expect_p2)
    a_equal(cube_gid[:,:,2], expect_p3)
    a_equal(cube_gid[:,:,3], expect_p4)
    a_equal(cube_gid[:,:,4], expect_p5)
    a_equal(cube_gid[:,:,5], expect_p6)



def test_make_global_sfc_3():
    '''
    cube_partition_sfc: make_global_sfc(ne=3)
    '''
    expect_p1 = np.array( \
            [[ 3, 2, 1], \
             [ 4, 7, 8], \
             [ 5, 6, 9]])

    expect_p2 = np.array( \
            [[12,11,10], \
             [13,16,17], \
             [14,15,18]])

    expect_p3 = np.array( \
            [[46,47,48], \
             [53,52,49], \
             [54,51,50]])

    expect_p4 = np.array( \
            [[36,35,28], \
             [33,34,29], \
             [32,31,30]])

    expect_p5 = np.array( \
            [[37,38,39], \
             [44,43,40], \
             [45,42,41]])

    expect_p6 = np.array( \
            [[23,24,27], \
             [22,25,26], \
             [21,20,19]])

    obj = CubePartitionSFC(ne=3, nproc=1)
    cube_gid = obj.make_global_sfc()
    a_equal(cube_gid[:,:,0], expect_p1)
    a_equal(cube_gid[:,:,1], expect_p2)
    a_equal(cube_gid[:,:,2], expect_p3)
    a_equal(cube_gid[:,:,3], expect_p4)
    a_equal(cube_gid[:,:,4], expect_p5)
    a_equal(cube_gid[:,:,5], expect_p6)



def test_make_cube_rank_2_1():
    '''
    cube_partition_sfc: make_cube_rank(ne=2, nproc=1)
    '''
    obj = CubePartitionSFC(ne=2, nproc=1)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(nelems, [24])
    a_equal(cube_rank, 0)
    a_equal(cube_lid[:,:,0], [[ 1, 3],
                              [ 2, 4]])
    a_equal(cube_lid[:,:,1], [[ 5, 7],
                              [ 6, 8]])
    a_equal(cube_lid[:,:,2], [[ 9,11],
                              [10,12]])
    a_equal(cube_lid[:,:,3], [[13,15],
                              [14,16]])
    a_equal(cube_lid[:,:,4], [[17,19],
                              [18,20]])
    a_equal(cube_lid[:,:,5], [[21,23],
                              [22,24]])



def test_make_cube_rank_2_8():
    '''
    cube_partition_sfc: make_cube_rank(ne=2, nproc=8)
    '''
    obj = CubePartitionSFC(ne=2, nproc=8)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(nelems, [3, 3, 3, 3, 3, 3, 3, 3])

    a_equal(cube_rank[:,:,0], [[ 0, 0],
                               [ 0, 1]])
    a_equal(cube_rank[:,:,1], [[ 1, 1],
                               [ 2, 2]])
    a_equal(cube_rank[:,:,2], [[ 6, 7],
                               [ 7, 7]])
    a_equal(cube_rank[:,:,3], [[ 5, 4],
                               [ 4, 4]])
    a_equal(cube_rank[:,:,4], [[ 5, 5],
                               [ 6, 6]])
    a_equal(cube_rank[:,:,5], [[ 3, 3],
                               [ 3, 2]])

    a_equal(cube_lid[:,:,0], [[ 1, 3],
                              [ 2, 1]])
    a_equal(cube_lid[:,:,1], [[ 2, 3],
                              [ 1, 2]])
    a_equal(cube_lid[:,:,2], [[ 1, 2],
                              [ 1, 3]])
    a_equal(cube_lid[:,:,3], [[ 1, 2],
                              [ 1, 3]])
    a_equal(cube_lid[:,:,4], [[ 2, 3],
                              [ 2, 3]])
    a_equal(cube_lid[:,:,5], [[ 1, 3],
                              [ 2, 3]])



def test_make_cube_rank_3_4():
    '''
    cube_partition_sfc: make_cube_rank(ne=3, nproc=4)
    '''
    obj = CubePartitionSFC(ne=3, nproc=4)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(nelems, [14, 14, 13, 13])

    a_equal(cube_rank[:,:,0], [[ 0, 0, 0],
                               [ 0, 0, 0],
                               [ 0, 0, 0]])
    a_equal(cube_rank[:,:,1], [[ 0, 0, 0],
                               [ 0, 1, 1],
                               [ 0, 1, 1]])
    a_equal(cube_rank[:,:,2], [[ 3, 3, 3],
                               [ 3, 3, 3],
                               [ 3, 3, 3]])
    a_equal(cube_rank[:,:,3], [[ 2, 2, 1],
                               [ 2, 2, 2],
                               [ 2, 2, 2]])
    a_equal(cube_rank[:,:,4], [[ 2, 2, 2],
                               [ 3, 3, 2],
                               [ 3, 3, 2]])
    a_equal(cube_rank[:,:,5], [[ 1, 1, 1],
                               [ 1, 1, 1],
                               [ 1, 1, 1]])

    a_equal(cube_lid[:,:,0], [[ 1, 4, 7],
                              [ 2, 5, 8],
                              [ 3, 6, 9]])
    a_equal(cube_lid[:,:,1], [[10,13,14],
                              [11, 1, 3],
                              [12, 2, 4]])
    a_equal(cube_lid[:,:,2], [[ 1, 4, 7],
                              [ 2, 5, 8],
                              [ 3, 6, 9]])
    a_equal(cube_lid[:,:,3], [[ 1, 4, 5],
                              [ 2, 5, 7],
                              [ 3, 6, 8]])
    a_equal(cube_lid[:,:,4], [[ 9,10,11],
                              [10,12,12],
                              [11,13,13]])
    a_equal(cube_lid[:,:,5], [[ 6, 9,12],
                              [ 7,10,13],
                              [ 8,11,14]])


def test_make_cube_rank_3_7():
    '''
    cube_partition_sfc: make_cube_rank(ne=3, nproc=7)
    '''
    obj = CubePartitionSFC(ne=3, nproc=7)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(nelems, [8, 8, 8, 8, 8, 7, 7])

    a_equal(cube_rank[:,:,0], [[ 0, 0, 0],
                               [ 0, 0, 0],
                               [ 0, 0, 1]])
    a_equal(cube_rank[:,:,1], [[ 1, 1, 1],
                               [ 1, 1, 2],
                               [ 1, 1, 2]])
    a_equal(cube_rank[:,:,2], [[ 5, 5, 6],
                               [ 6, 6, 6],
                               [ 6, 6, 6]])
    a_equal(cube_rank[:,:,3], [[ 4, 4, 3],
                               [ 4, 4, 3],
                               [ 3, 3, 3]])
    a_equal(cube_rank[:,:,4], [[ 4, 4, 4],
                               [ 5, 5, 4],
                               [ 5, 5, 5]])
    a_equal(cube_rank[:,:,5], [[ 2, 2, 3],
                               [ 2, 3, 3],
                               [ 2, 2, 2]])

    a_equal(cube_lid[:,:,0], [[ 1, 4, 7],
                              [ 2, 5, 8],
                              [ 3, 6, 1]])
    a_equal(cube_lid[:,:,1], [[ 2, 5, 8],
                              [ 3, 6, 1],
                              [ 4, 7, 2]])
    a_equal(cube_lid[:,:,2], [[ 1, 2, 5],
                              [ 1, 3, 6],
                              [ 2, 4, 7]])
    a_equal(cube_lid[:,:,3], [[ 1, 3, 3],
                              [ 2, 4, 4],
                              [ 1, 2, 5]])
    a_equal(cube_lid[:,:,4], [[ 5, 6, 7],
                              [ 3, 5, 8],
                              [ 4, 6, 7]])
    a_equal(cube_lid[:,:,5], [[ 3, 6, 7],
                              [ 4, 6, 8],
                              [ 5, 7, 8]])



def test_make_cube_rank_4_5():
    '''
    cube_partition_sfc: make_cube_rank(ne=4, nproc=5)
    '''
    obj = CubePartitionSFC(ne=4, nproc=5)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(nelems, [20, 19, 19, 19, 19])

    a_equal(cube_rank[:,:,0], [[ 0, 0, 0, 0],
                               [ 0, 0, 0, 0],
                               [ 0, 0, 0, 0],
                               [ 0, 0, 0, 0]])
    a_equal(cube_rank[:,:,1], [[ 1, 1, 0, 0],
                               [ 1, 1, 0, 0],
                               [ 1, 1, 1, 1],
                               [ 1, 1, 1, 1]])
    a_equal(cube_rank[:,:,2], [[ 4, 4, 4, 4],
                               [ 4, 4, 4, 4],
                               [ 4, 4, 4, 4],
                               [ 4, 4, 4, 4]])
    a_equal(cube_rank[:,:,3], [[ 3, 3, 2, 2],
                               [ 3, 3, 2, 2],
                               [ 3, 2, 2, 2],
                               [ 3, 2, 2, 2]])
    a_equal(cube_rank[:,:,4], [[ 3, 3, 3, 3],
                               [ 3, 3, 3, 3],
                               [ 4, 4, 3, 3],
                               [ 4, 3, 3, 3]])
    a_equal(cube_rank[:,:,5], [[ 2, 2, 2, 2],
                               [ 2, 2, 2, 2],
                               [ 1, 2, 1, 1],
                               [ 1, 1, 1, 1]])

    a_equal(cube_lid[:,:,0], [[ 1, 5, 9,13],
                              [ 2, 6,10,14],
                              [ 3, 7,11,15],
                              [ 4, 8,12,16]])
    a_equal(cube_lid[:,:,1], [[ 1, 5,17,19],
                              [ 2, 6,18,20],
                              [ 3, 7, 9,11],
                              [ 4, 8,10,12]])
    a_equal(cube_lid[:,:,2], [[ 1, 5, 9,13],
                              [ 2, 6,10,14],
                              [ 3, 7,11,15],
                              [ 4, 8,12,16]])
    a_equal(cube_lid[:,:,3], [[ 1, 5, 3, 7],
                              [ 2, 6, 4, 8],
                              [ 3, 1, 5, 9],
                              [ 4, 2, 6,10]])
    a_equal(cube_lid[:,:,4], [[ 7, 9,12,16],
                              [ 8,10,13,17],
                              [17,19,14,18],
                              [18,11,15,19]])
    a_equal(cube_lid[:,:,5], [[11,13,16,18],
                              [12,14,17,19],
                              [13,15,16,18],
                              [14,15,17,19]])



def test_make_cube_rank_6_10():
    '''
    cube_partition_sfc: make_cube_rank(ne=6, nproc=10)
    '''
    obj = CubePartitionSFC(ne=6, nproc=10)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(cube_rank[:,:,0], [[ 0, 0, 0, 0, 0, 0],
                               [ 0, 0, 0, 0, 0, 0],
                               [ 0, 0, 1, 1, 1, 1],
                               [ 0, 0, 1, 1, 1, 1],
                               [ 0, 0, 1, 1, 1, 1],
                               [ 0, 0, 0, 0, 1, 1]])
    a_equal(cube_rank[:,:,1], [[ 2, 2, 1, 1, 1, 1],
                               [ 2, 2, 1, 1, 1, 1],
                               [ 2, 2, 2, 2, 2, 3],
                               [ 2, 2, 2, 2, 2, 3],
                               [ 2, 2, 2, 2, 3, 3],
                               [ 2, 2, 2, 2, 3, 3]])
    a_equal(cube_rank[:,:,2], [[ 8, 8, 8, 8, 8, 8],
                               [ 8, 8, 8, 8, 8, 8],
                               [ 9, 9, 9, 9, 8, 8],
                               [ 9, 9, 9, 9, 9, 8],
                               [ 9, 9, 9, 9, 9, 9],
                               [ 9, 9, 9, 9, 9, 9]])
    a_equal(cube_rank[:,:,3], [[ 6, 6, 6, 6, 4, 4],
                               [ 6, 6, 6, 6, 5, 5],
                               [ 5, 5, 6, 6, 5, 5],
                               [ 5, 5, 6, 6, 5, 5],
                               [ 5, 5, 5, 5, 5, 5],
                               [ 5, 5, 5, 5, 5, 5]])
    a_equal(cube_rank[:,:,4], [[ 6, 6, 6, 6, 6, 7],
                               [ 6, 6, 6, 6, 7, 7],
                               [ 8, 7, 7, 7, 7, 7],
                               [ 8, 7, 7, 7, 7, 7],
                               [ 8, 8, 7, 7, 7, 7],
                               [ 8, 8, 7, 7, 7, 7]])
    a_equal(cube_rank[:,:,5], [[ 4, 4, 4, 4, 4, 4],
                               [ 4, 4, 4, 4, 4, 4],
                               [ 3, 3, 4, 4, 4, 4],
                               [ 3, 3, 4, 4, 4, 4],
                               [ 3, 3, 3, 3, 3, 3],
                               [ 3, 3, 3, 3, 3, 3]])

    a_equal(cube_lid[:,:,0], [[ 1, 7,13,16,19,21],
                              [ 2, 8,14,17,20,22],
                              [ 3, 9, 1, 4, 7,11],
                              [ 4,10, 2, 5, 8,12],
                              [ 5,11, 3, 6, 9,13],
                              [ 6,12,15,18,10,14]])
    a_equal(cube_lid[:,:,1], [[ 1, 7,15,17,19,21],
                              [ 2, 8,16,18,20,22],
                              [ 3, 9,13,17,21, 3],
                              [ 4,10,14,18,22, 4],
                              [ 5,11,15,19, 1, 5],
                              [ 6,12,16,20, 2, 6]])
    a_equal(cube_lid[:,:,2], [[ 1, 3, 5, 7, 9,12],
                              [ 2, 4, 6, 8,10,13],
                              [ 1, 5, 9,13,11,14],
                              [ 2, 6,10,14,17,15],
                              [ 3, 7,11,15,18,20],
                              [ 4, 8,12,16,19,21]])
    a_equal(cube_lid[:,:,3], [[ 1, 3, 5, 9, 1, 2],
                              [ 2, 4, 6,10,13,18],
                              [ 1, 5, 7,11,14,19],
                              [ 2, 6, 8,12,15,20],
                              [ 3, 7, 9,11,16,21],
                              [ 4, 8,10,12,17,22]])
    a_equal(cube_lid[:,:,4], [[13,15,17,19,21,16],
                              [14,16,18,20,11,17],
                              [16, 1, 3, 7,12,18],
                              [17, 2, 4, 8,13,19],
                              [18,20, 5, 9,14,20],
                              [19,21, 6,10,15,21]])
    a_equal(cube_lid[:,:,5], [[ 3, 5, 7,11,15,19],
                              [ 4, 6, 8,12,16,20],
                              [ 7,11, 9,13,17,21],
                              [ 8,12,10,14,18,22],
                              [ 9,13,15,17,19,21],
                              [10,14,16,18,20,22]])


"""
def test_make_cube_rank_30_24():
    '''
    cube_partition_sfc: make_cube_rank(ne=30, nproc=24)
    '''
    obj = CubePartitionSFC(ne=30, nproc=24)
    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    a_equal(nelems, 225)

    a_equal(cube_rank[:15,15:,0], 0)
    a_equal(cube_rank[:15,:15,0], 1)
    a_equal(cube_rank[15:,:15,0], 2)
    a_equal(cube_rank[15:,15:,0], 3)

    a_equal(cube_rank[:15,15:,1], 4)
    a_equal(cube_rank[:15,:15,1], 5)
    a_equal(cube_rank[15:,:15,1], 6)
    a_equal(cube_rank[15:,15:,1], 7)

    a_equal(cube_rank[15:,15:,5], 8)
    a_equal(cube_rank[15:,:15,5], 9)
    a_equal(cube_rank[:15,:15,5],10)
    a_equal(cube_rank[:15,15:,5],11)

    a_equal(cube_rank[:15,15:,3],12)
    a_equal(cube_rank[15:,15:,3],13)
    a_equal(cube_rank[15:,:15,3],14)
    a_equal(cube_rank[:15,:15,3],15)

    a_equal(cube_rank[:15,:15,4],16)
    a_equal(cube_rank[:15,15:,4],17)
    a_equal(cube_rank[15:,15:,4],18)
    a_equal(cube_rank[15:,:15,4],19)

    a_equal(cube_rank[:15,:15,2],20)
    a_equal(cube_rank[:15,15:,2],21)
    a_equal(cube_rank[15:,15:,2],22)
    a_equal(cube_rank[15:,:15,2],23)
"""



def test_make_elem_coord_6_10():
    '''
    cube_partition_sfc: make_elem_coord(): ne=6, nproc=10
    '''
    obj = CubePartitionSFC(ne=6, nproc=10)

    nelems, cube_rank, cube_lid = obj.make_cube_rank()
    iproc = 0
    nelem = nelems[0]
    elem_coord = obj.make_elem_coord(iproc, nelem, cube_rank, cube_lid)
    a_equal(elem_coord, np.array(
        [ 1, 1, 1,
          2, 1, 1,
          3, 1, 1,
          4, 1, 1,
          5, 1, 1,
          6, 1, 1,
          1, 2, 1,
          2, 2, 1,
          3, 2, 1,
          4, 2, 1,
          5, 2, 1,
          6, 2, 1,
          1, 3, 1,
          2, 3, 1,
          6, 3, 1,
          1, 4, 1,
          2, 4, 1,
          6, 4, 1,
          1, 5, 1,
          2, 5, 1,
          1, 6, 1,
          2, 6, 1], 'i4').reshape((3,22), order='F'))

    elem_coord = obj.make_elem_coord(1, nelems[1], cube_rank, cube_lid)
    a_equal(elem_coord, np.array(
        [ 3, 3, 1,
          4, 3, 1,
          5, 3, 1,
          3, 4, 1,
          4, 4, 1,
          5, 4, 1,
          3, 5, 1,
          4, 5, 1,
          5, 5, 1,
          6, 5, 1,
          3, 6, 1,
          4, 6, 1,
          5, 6, 1,
          6, 6, 1,
          1, 3, 2,
          2, 3, 2,
          1, 4, 2,
          2, 4, 2,
          1, 5, 2,
          2, 5, 2,
          1, 6, 2,
          2, 6, 2], 'i4').reshape((3,22), order='F'))

    elem_coord = obj.make_elem_coord(9, nelems[9], cube_rank, cube_lid)
    a_equal(elem_coord, np.array(
        [ 3, 1, 3,
          4, 1, 3,
          5, 1, 3,
          6, 1, 3,
          3, 2, 3,
          4, 2, 3,
          5, 2, 3,
          6, 2, 3,
          3, 3, 3,
          4, 3, 3,
          5, 3, 3,
          6, 3, 3,
          3, 4, 3,
          4, 4, 3,
          5, 4, 3,
          6, 4, 3,
          4, 5, 3,
          5, 5, 3,
          6, 5, 3,
          5, 6, 3,
          6, 6, 3], 'i4').reshape((3,21), order='F'))



#==============================================================================
#==============================================================================



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



def plot_cube_partition(ne, nproc, save):
    '''
    cube_partition_sfc: plot the cube_rank array
    '''
    obj = CubePartitionSFC(ne, nproc)

    nelems, cube_rank, cube_lid = obj.make_cube_rank()

    obj2 = CubePartitionBand(ne, nproc)
    cube_color = obj2.make_cube_color(cube_rank)
    perimeter_ratio, num_nbrs = obj2.global_perimeter_ratio(cube_rank)
    print('ne={}, nproc={}, perimter_ratio={}'.format(ne, nproc, perimeter_ratio))

    box = np.ones((4*ne,3*ne), 'i4', order='F')*(-1)
    box[:ne,2*ne:] = np.rot90(np.fliplr(cube_color[:,:,0]), 3)
    box[:ne,ne:2*ne] = np.rot90(np.fliplr(cube_color[:,:,1]), 3)
    box[ne:2*ne,ne:2*ne] = cube_color[:,:,5]
    box[2*ne:3*ne,ne:2*ne] = np.rot90(cube_color[:,:,3], 2)
    box[3*ne:,ne:2*ne] = cube_color[:,:,4]
    box[3*ne:,:ne] = cube_color[:,:,2]

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
    boxes = [mp.Rectangle((-0.5, 2*ne-0.5), ne, ne, **kwds),
             mp.Rectangle((-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((ne-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((2*ne-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((3*ne-0.5, ne-0.5), ne, ne, **kwds),
             mp.Rectangle((3*ne-0.5, -0.5), ne, ne, **kwds)]
    for box in boxes: ax.add_patch(box)

    plt.tight_layout()
    if save:
        fname = 'cube_partition_sfc.ne{}_nproc{}.png'.format(ne, nproc)
        plt.savefig('png/'+fname, dpi=300)
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--save', action='store_true', help='save as png format')
    parser.add_argument('ne', type=int, help='number of elements')
    parser.add_argument('nproc', type=int, help='number of processors')
    args = parser.parse_args()

    plot_cube_partition(args.ne, args.nproc, args.save)
