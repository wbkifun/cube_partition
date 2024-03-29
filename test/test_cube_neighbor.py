'''

abstract : unittest of cube_neighbor.f90

history :
  2017-04-15  ki-hwan kim  start
  2017-05-09  ki-hwan kim  remove convert_neighbor_ij(), find_nbr_elem_ij()
                           add convert_nbr_elem_ij()
  2017-05-15  ki-hwan kim  convert_nbr_elem_ij() -> convert_nbr_eij()
                           add convert_nbr_gij()

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join
import sys

from numpy.testing import assert_equal as equal
from numpy.testing import assert_array_equal as a_equal
import numpy as np
import matplotlib.pyplot as plt


current_dir = dirname(abspath(__file__))
sys.path.append(dirname(current_dir))
from cube_neighbor import CubeNeighbor



obj = CubeNeighbor()



def test_quotient():
    '''
    cube_neighbor: quotient(): n=3
    '''
    n = 3
    i_list = list(range(-5,8))
    i_list = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
    q_list = [-2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2]
    for i, q in zip(i_list, q_list):
        equal(obj.quotient(n, i), q)


def test_convert_rotated_ij_6_2_3():
    '''
    cube_neighbor: convert_rotated_ij(ne=6, ei=2, ej=3)
    '''
    ne = 6
    ei, ej = 2, 3
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 0), (2, 3))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 1), (3, 5))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 2), (5, 4))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 3), (4, 2))



def test_convert_rotated_ij_6_1_4():
    '''
    cube_neighbor: convert_rotated_ij(ne=6, ei=1, ej=4)
    '''
    ne = 6
    ei, ej = 1, 4
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 0), (1, 4))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 1), (4, 6))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 2), (6, 3))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 3), (3, 1))



def test_convert_rotated_ij_6_5_1():
    '''
    cube_neighbor: convert_rotated_ij(ne=6, ei=5, ej=1)
    '''
    ne = 6
    ei, ej = 5, 1
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 0), (5, 1))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 1), (1, 2))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 2), (2, 6))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 3), (6, 5))



def test_convert_rotated_ij_6_0_4():
    '''
    cube_neighbor: convert_rotated_ij(ne=6, ei=0, ej=4), out of bound
    '''
    ne = 6
    ei, ej = 0, 4
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 0), (0, 4))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 1), (4, 7))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 2), (7, 3))
    a_equal(obj.convert_rotated_ij(ne, ei, ej, 3), (3, 0))



def test_convert_nbr_eij_6_1():
    '''
    cube_neighbor: convert_nbr_eij(): ne=6, panel=1
    '''
    ne = 6

    a_equal(obj.convert_nbr_eij(ne, ei=1, ej=3, panel=1), (1, 3, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, 6, 3, 1), (6, 3, 1, 0))

    a_equal(obj.convert_nbr_eij(ne,  7, 3, 1), (1, 3, 2, 0))
    a_equal(obj.convert_nbr_eij(ne,  8, 3, 1), (2, 3, 2, 0))
    a_equal(obj.convert_nbr_eij(ne, 12, 3, 1), (6, 3, 2, 0))
    #a_equal(obj.convert_nbr_eij(ne, 13, 3, 1), (1, 3, 3, 0))

    a_equal(obj.convert_nbr_eij(ne,  0, 3, 1), (6, 3, 4, 0))
    a_equal(obj.convert_nbr_eij(ne, -1, 3, 1), (5, 3, 4, 0))
    a_equal(obj.convert_nbr_eij(ne, -5, 3, 1), (1, 3, 4, 0))
    #a_equal(obj.convert_nbr_eij(ne, -6, 3, 1), (6, 3, 3, 0))

    a_equal(obj.convert_nbr_eij(ne, 1, 7, 1), (1, 1, 6, 0))
    a_equal(obj.convert_nbr_eij(ne, 1, 8, 1), (1, 2, 6, 0))
    a_equal(obj.convert_nbr_eij(ne, 1,12, 1), (1, 6, 6, 0))
    #a_equal(obj.convert_nbr_eij(ne, 1,13, 1), (3, 6, 3, 2))

    a_equal(obj.convert_nbr_eij(ne, 1, 0, 1), (1, 6, 5, 0))
    a_equal(obj.convert_nbr_eij(ne, 1,-1, 1), (1, 5, 5, 0))
    a_equal(obj.convert_nbr_eij(ne, 1,-5, 1), (1, 1, 5, 0))
    #a_equal(obj.convert_nbr_eij(ne, 1,-6, 1), (3, 1, 3, 2))
    
    a_equal(obj.convert_nbr_eij(ne, 0, 0, 1), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 0, 1), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 0, 7, 1), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 7, 1), (-1, -1, -1, -1))



def test_convert_nbr_eij_6_2():
    '''
    cube_neighbor: convert_nbr_eij(): ne=6, panel=2
    '''
    ne = 6

    a_equal(obj.convert_nbr_eij(ne, ei=1, ej=3, panel=2), (1, 3, 2, 0))
    a_equal(obj.convert_nbr_eij(ne, 6, 3, 2), (6, 3, 2, 0))

    a_equal(obj.convert_nbr_eij(ne,  7, 3, 2), (1, 3, 3, 0))
    a_equal(obj.convert_nbr_eij(ne,  8, 3, 2), (2, 3, 3, 0))
    a_equal(obj.convert_nbr_eij(ne, 12, 3, 2), (6, 3, 3, 0))
    #a_equal(obj.convert_nbr_eij(ne, 13, 3, 2), (1, 3, 4, 0))

    a_equal(obj.convert_nbr_eij(ne,  0, 3, 2), (6, 3, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, -1, 3, 2), (5, 3, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, -5, 3, 2), (1, 3, 1, 0))
    #a_equal(obj.convert_nbr_eij(ne, -6, 3, 2), (6, 3, 4, 0))

    a_equal(obj.convert_nbr_eij(ne, 1, 7, 2), (6, 1, 6, 3))
    a_equal(obj.convert_nbr_eij(ne, 1, 8, 2), (5, 1, 6, 3))
    a_equal(obj.convert_nbr_eij(ne, 1,12, 2), (1, 1, 6, 3))
    #a_equal(obj.convert_nbr_eij(ne, 1,13, 2), (6, 6, 4, 2))

    a_equal(obj.convert_nbr_eij(ne, 1, 0, 2), (6, 6, 5, 1))
    a_equal(obj.convert_nbr_eij(ne, 1,-1, 2), (5, 6, 5, 1))
    a_equal(obj.convert_nbr_eij(ne, 1,-5, 2), (1, 6, 5, 1))
    #a_equal(obj.convert_nbr_eij(ne, 1,-6, 2), (3, 1, 4, 2))
    
    a_equal(obj.convert_nbr_eij(ne, 0, 0, 2), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 0, 2), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 0, 7, 2), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 7, 2), (-1, -1, -1, -1))



def test_convert_nbr_eij_6_3():
    '''
    cube_neighbor: convert_nbr_eij(): ne=6, panel=3
    '''
    ne = 6

    a_equal(obj.convert_nbr_eij(ne, ei=1, ej=3, panel=3), (1, 3, 3, 0))
    a_equal(obj.convert_nbr_eij(ne, 6, 3, 3), (6, 3, 3, 0))

    a_equal(obj.convert_nbr_eij(ne,  7, 3, 3), (1, 3, 4, 0))
    a_equal(obj.convert_nbr_eij(ne,  8, 3, 3), (2, 3, 4, 0))
    a_equal(obj.convert_nbr_eij(ne, 12, 3, 3), (6, 3, 4, 0))
    #a_equal(obj.convert_nbr_eij(ne, 13, 3, 3), (1, 3, 1, 0))

    a_equal(obj.convert_nbr_eij(ne,  0, 3, 3), (6, 3, 2, 0))
    a_equal(obj.convert_nbr_eij(ne, -1, 3, 3), (5, 3, 2, 0))
    a_equal(obj.convert_nbr_eij(ne, -5, 3, 3), (1, 3, 2, 0))
    #a_equal(obj.convert_nbr_eij(ne, -6, 3, 3), (6, 3, 1, 0))

    a_equal(obj.convert_nbr_eij(ne, 1, 7, 3), (6, 6, 6, 2))
    a_equal(obj.convert_nbr_eij(ne, 1, 8, 3), (6, 5, 6, 2))
    a_equal(obj.convert_nbr_eij(ne, 1,12, 3), (6, 1, 6, 2))
    #a_equal(obj.convert_nbr_eij(ne, 1,13, 3), (6, 6, 1, 2))

    a_equal(obj.convert_nbr_eij(ne, 1, 0, 3), (6, 1, 5, 2))
    a_equal(obj.convert_nbr_eij(ne, 1,-1, 3), (6, 2, 5, 2))
    a_equal(obj.convert_nbr_eij(ne, 1,-5, 3), (6, 6, 5, 2))
    #a_equal(obj.convert_nbr_eij(ne, 1,-6, 3), (6, 1, 1, 2))
    
    a_equal(obj.convert_nbr_eij(ne, 0, 0, 3), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 0, 3), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 0, 7, 3), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 7, 3), (-1, -1, -1, -1))



def test_convert_nbr_eij_6_4():
    '''
    cube_neighbor: convert_nbr_eij(): ne=6, panel=4
    '''
    ne = 6

    a_equal(obj.convert_nbr_eij(ne, ei=1, ej=3, panel=4), (1, 3, 4, 0))
    a_equal(obj.convert_nbr_eij(ne, 6, 3, 4), (6, 3, 4, 0))

    a_equal(obj.convert_nbr_eij(ne,  7, 3, 4), (1, 3, 1, 0))
    a_equal(obj.convert_nbr_eij(ne,  8, 3, 4), (2, 3, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, 12, 3, 4), (6, 3, 1, 0))
    #a_equal(obj.convert_nbr_eij(ne, 13, 3, 4), (1, 3, 2, 0))

    a_equal(obj.convert_nbr_eij(ne,  0, 3, 4), (6, 3, 3, 0))
    a_equal(obj.convert_nbr_eij(ne, -1, 3, 4), (5, 3, 3, 0))
    a_equal(obj.convert_nbr_eij(ne, -5, 3, 4), (1, 3, 3, 0))
    #a_equal(obj.convert_nbr_eij(ne, -6, 3, 4), (6, 3, 2, 0))

    a_equal(obj.convert_nbr_eij(ne, 1, 7, 4), (1, 6, 6, 1))
    a_equal(obj.convert_nbr_eij(ne, 1, 8, 4), (2, 6, 6, 1))
    a_equal(obj.convert_nbr_eij(ne, 1,12, 4), (6, 6, 6, 1))
    #a_equal(obj.convert_nbr_eij(ne, 1,13, 4), (6, 6, 2, 2))

    a_equal(obj.convert_nbr_eij(ne, 1, 0, 4), (1, 1, 5, 3))
    a_equal(obj.convert_nbr_eij(ne, 1,-1, 4), (2, 1, 5, 3))
    a_equal(obj.convert_nbr_eij(ne, 1,-5, 4), (6, 1, 5, 3))
    #a_equal(obj.convert_nbr_eij(ne, 1,-6, 4), (6, 1, 2, 2))
    
    a_equal(obj.convert_nbr_eij(ne, 0, 0, 4), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 0, 4), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 0, 7, 4), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 7, 4), (-1, -1, -1, -1))



def test_convert_nbr_eij_6_5():
    '''
    cube_neighbor: convert_nbr_eij(): ne=6, panel=5
    '''
    ne = 6

    a_equal(obj.convert_nbr_eij(ne, ei=1, ej=3, panel=5), (1, 3, 5, 0))
    a_equal(obj.convert_nbr_eij(ne, 6, 3, 5), (6, 3, 5, 0))

    a_equal(obj.convert_nbr_eij(ne,  7, 3, 5), (4, 1, 2, 3))
    a_equal(obj.convert_nbr_eij(ne,  8, 3, 5), (4, 2, 2, 3))
    a_equal(obj.convert_nbr_eij(ne, 12, 3, 5), (4, 6, 2, 3))
    #a_equal(obj.convert_nbr_eij(ne, 13, 3, 5), (6, 4, 6, 2))

    a_equal(obj.convert_nbr_eij(ne,  0, 3, 5), (3, 1, 4, 1))
    a_equal(obj.convert_nbr_eij(ne, -1, 3, 5), (3, 2, 4, 1))
    a_equal(obj.convert_nbr_eij(ne, -5, 3, 5), (3, 6, 4, 1))
    #a_equal(obj.convert_nbr_eij(ne, -6, 3, 5), (1, 4, 6, 2))

    a_equal(obj.convert_nbr_eij(ne, 1, 7, 5), (1, 1, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, 1, 8, 5), (1, 2, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, 1,12, 5), (1, 6, 1, 0))
    #a_equal(obj.convert_nbr_eij(ne, 1,13, 5), (1, 1, 6, 0))

    a_equal(obj.convert_nbr_eij(ne, 1, 0, 5), (6, 1, 3, 2))
    a_equal(obj.convert_nbr_eij(ne, 1,-1, 5), (6, 2, 3, 2))
    a_equal(obj.convert_nbr_eij(ne, 1,-5, 5), (6, 6, 3, 2))
    #a_equal(obj.convert_nbr_eij(ne, 1,-6, 5), (1, 6, 6, 0))
    
    a_equal(obj.convert_nbr_eij(ne, 0, 0, 5), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 0, 5), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 0, 7, 5), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 7, 5), (-1, -1, -1, -1))

    # for paper
    ne = 3
    a_equal(obj.convert_nbr_eij(ne, 3, -1, 5), (1, 2, 3, 2))



def test_convert_nbr_eij_6_6():
    '''
    cube_neighbor: convert_nbr_eij(): ne=6, panel=6
    '''
    ne = 6

    a_equal(obj.convert_nbr_eij(ne, ei=1, ej=3, panel=6), (1, 3, 6, 0))
    a_equal(obj.convert_nbr_eij(ne, 6, 3, 6), (6, 3, 6, 0))

    a_equal(obj.convert_nbr_eij(ne,  7, 3, 6), (3, 6, 2, 1))
    a_equal(obj.convert_nbr_eij(ne,  8, 3, 6), (3, 5, 2, 1))
    a_equal(obj.convert_nbr_eij(ne, 12, 3, 6), (3, 1, 2, 1))
    #a_equal(obj.convert_nbr_eij(ne, 13, 3, 6), (6, 4, 5, 2))

    a_equal(obj.convert_nbr_eij(ne,  0, 3, 6), (4, 6, 4, 3))
    a_equal(obj.convert_nbr_eij(ne, -1, 3, 6), (4, 5, 4, 3))
    a_equal(obj.convert_nbr_eij(ne, -5, 3, 6), (4, 1, 4, 3))
    #a_equal(obj.convert_nbr_eij(ne, -6, 3, 6), (1, 4, 5, 2))

    a_equal(obj.convert_nbr_eij(ne, 1, 7, 6), (6, 6, 3, 2))
    a_equal(obj.convert_nbr_eij(ne, 1, 8, 6), (6, 5, 3, 2))
    a_equal(obj.convert_nbr_eij(ne, 1,12, 6), (6, 1, 3, 2))
    #a_equal(obj.convert_nbr_eij(ne, 1,13, 6), (1, 1, 5, 0))

    a_equal(obj.convert_nbr_eij(ne, 1, 0, 6), (1, 6, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, 1,-1, 6), (1, 5, 1, 0))
    a_equal(obj.convert_nbr_eij(ne, 1,-5, 6), (1, 1, 1, 0))
    #a_equal(obj.convert_nbr_eij(ne, 1,-6, 6), (1, 6, 5, 0))
    
    a_equal(obj.convert_nbr_eij(ne, 0, 0, 6), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 0, 6), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 0, 7, 6), (-1, -1, -1, -1))
    a_equal(obj.convert_nbr_eij(ne, 7, 7, 6), (-1, -1, -1, -1))



def test_convert_nbr_eij_10_4():
    '''
    cube_neighbor: convert_nbr_eij(): ne=10, panel=4
    '''
    ne = 10
    a_equal(obj.convert_nbr_eij(ne, ei=6, ej=10, panel=4), (6, 10, 4, 0))
    a_equal(obj.convert_nbr_eij(ne, ei=6, ej=11, panel=4), (1, 5, 6, 1))



def test_convert_nbr_gij_6_1():
    '''
    cube_neighbor: convert_nbr_gij(): ne=6, np=4, panel=1
    '''
    ne = 6
    ngq = 4

    a_equal(obj.convert_nbr_gij(ne, ngq, gi=1, gj=1, ei=1, ej=3, panel=1), (1, 1, 1, 3, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 4, 4, 6, 3, 1), (4, 4, 6, 3, 1))

    a_equal(obj.convert_nbr_gij(ne, ngq,  8, 1, 6, 3, 1), (4, 1, 1, 3, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 12, 1, 6, 3, 1), (4, 1, 2, 3, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 28, 1, 6, 3, 1), (4, 1, 6, 3, 2))

    a_equal(obj.convert_nbr_gij(ne, ngq,   0, 1, 1, 3, 1), (4, 1, 6, 3, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq,  -4, 1, 1, 3, 1), (4, 1, 5, 3, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq, -20, 1, 1, 3, 1), (4, 1, 1, 3, 4))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1, 8, 1, 6, 1), (1, 4, 1, 1, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,12, 1, 6, 1), (1, 4, 1, 2, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,28, 1, 6, 1), (1, 4, 1, 6, 6))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1,  0, 1, 1, 1), (1, 4, 1, 6, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1, -4, 1, 1, 1), (1, 4, 1, 5, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,-20, 1, 1, 1), (1, 4, 1, 1, 5))
    
    a_equal(obj.convert_nbr_gij(ne, ngq,  0, 0, 1, 1, 1), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25, 0, 1, 1, 1), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  0,25, 1, 1, 1), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25,25, 1, 1, 1), (-1, -1, -1, -1, -1))



def test_convert_nbr_gij_6_2():
    '''
    cube_neighbor: convert_nbr_gij(): ne=6, np=4, panel=2
    '''
    ne = 6
    ngq = 4

    a_equal(obj.convert_nbr_gij(ne, ngq, gi=1, gj=1, ei=1, ej=3, panel=2), (1, 1, 1, 3, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 4, 4, 6, 3, 2), (4, 4, 6, 3, 2))

    a_equal(obj.convert_nbr_gij(ne, ngq,  8, 1, 6, 3, 2), (4, 1, 1, 3, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 12, 1, 6, 3, 2), (4, 1, 2, 3, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 28, 1, 6, 3, 2), (4, 1, 6, 3, 3))

    a_equal(obj.convert_nbr_gij(ne, ngq,   0, 1, 1, 3, 2), (4, 1, 6, 3, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  -4, 1, 1, 3, 2), (4, 1, 5, 3, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, -20, 1, 1, 3, 2), (4, 1, 1, 3, 1))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1, 8, 1, 6, 2), (1, 1, 6, 1, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,12, 1, 6, 2), (1, 1, 5, 1, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,28, 1, 6, 2), (1, 1, 1, 1, 6))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1,  0, 1, 1, 2), (4, 4, 6, 6, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1, -4, 1, 1, 2), (4, 4, 5, 6, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,-20, 1, 1, 2), (4, 4, 1, 6, 5))
    
    a_equal(obj.convert_nbr_gij(ne, ngq,  0, 0, 1, 1, 2), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25, 0, 1, 1, 2), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  0,25, 1, 1, 2), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25,25, 1, 1, 2), (-1, -1, -1, -1, -1))



def test_convert_nbr_gij_6_3():
    '''
    cube_neighbor: convert_nbr_gij(): ne=6, np=4, panel=3
    '''
    ne = 6
    ngq = 4

    a_equal(obj.convert_nbr_gij(ne, ngq, gi=1, gj=1, ei=1, ej=3, panel=3), (1, 1, 1, 3, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 4, 4, 6, 3, 3), (4, 4, 6, 3, 3))

    a_equal(obj.convert_nbr_gij(ne, ngq,  8, 1, 6, 3, 3), (4, 1, 1, 3, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq, 12, 1, 6, 3, 3), (4, 1, 2, 3, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq, 28, 1, 6, 3, 3), (4, 1, 6, 3, 4))

    a_equal(obj.convert_nbr_gij(ne, ngq,   0, 1, 1, 3, 3), (4, 1, 6, 3, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq,  -4, 1, 1, 3, 3), (4, 1, 5, 3, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, -20, 1, 1, 3, 3), (4, 1, 1, 3, 2))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1, 8, 1, 6, 3), (4, 1, 6, 6, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,12, 1, 6, 3), (4, 1, 6, 5, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,28, 1, 6, 3), (4, 1, 6, 1, 6))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1,  0, 1, 1, 3), (4, 1, 6, 1, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1, -4, 1, 1, 3), (4, 1, 6, 2, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,-20, 1, 1, 3), (4, 1, 6, 6, 5))
    
    a_equal(obj.convert_nbr_gij(ne, ngq,  0, 0, 1, 1, 3), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25, 0, 1, 1, 3), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  0,25, 1, 1, 3), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25,25, 1, 1, 3), (-1, -1, -1, -1, -1))



def test_convert_nbr_eij_6_4():
    '''
    cube_neighbor: convert_nbr_gij(): ne=6, np=4, panel=4
    '''
    ne = 6
    ngq = 4

    a_equal(obj.convert_nbr_gij(ne, ngq, gi=1, gj=1, ei=1, ej=3, panel=4), (1, 1, 1, 3, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq, 4, 4, 6, 3, 4), (4, 4, 6, 3, 4))

    a_equal(obj.convert_nbr_gij(ne, ngq,  8, 1, 6, 3, 4), (4, 1, 1, 3, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 12, 1, 6, 3, 4), (4, 1, 2, 3, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 28, 1, 6, 3, 4), (4, 1, 6, 3, 1))

    a_equal(obj.convert_nbr_gij(ne, ngq,   0, 1, 1, 3, 4), (4, 1, 6, 3, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq,  -4, 1, 1, 3, 4), (4, 1, 5, 3, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, -20, 1, 1, 3, 4), (4, 1, 1, 3, 3))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1, 8, 1, 6, 4), (4, 4, 1, 6, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,12, 1, 6, 4), (4, 4, 2, 6, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,28, 1, 6, 4), (4, 4, 6, 6, 6))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1,  0, 1, 1, 4), (1, 1, 1, 1, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1, -4, 1, 1, 4), (1, 1, 2, 1, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,-20, 1, 1, 4), (1, 1, 6, 1, 5))
    
    a_equal(obj.convert_nbr_gij(ne, ngq,  0, 0, 1, 1, 4), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25, 0, 1, 1, 4), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  0,25, 1, 1, 4), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25,25, 1, 1, 4), (-1, -1, -1, -1, -1))



def test_convert_nbr_gij_6_5():
    '''
    cube_neighbor: convert_nbr_gij(): ne=6, np=4, panel=5
    '''
    ne = 6
    ngq = 4

    a_equal(obj.convert_nbr_gij(ne, ngq, gi=1, gj=1, ei=1, ej=3, panel=5), (1, 1, 1, 3, 5))
    a_equal(obj.convert_nbr_gij(ne, ngq, 4, 4, 6, 3, 5), (4, 4, 6, 3, 5))

    a_equal(obj.convert_nbr_gij(ne, ngq,  8, 1, 6, 3, 5), (4, 4, 4, 1, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 12, 1, 6, 3, 5), (4, 4, 4, 2, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 28, 1, 6, 3, 5), (4, 4, 4, 6, 2))

    a_equal(obj.convert_nbr_gij(ne, ngq,   0, 1, 1, 3, 5), (1, 1, 3, 1, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq,  -4, 1, 1, 3, 5), (1, 1, 3, 2, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq, -20, 1, 1, 3, 5), (1, 1, 3, 6, 4))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1, 8, 1, 6, 5), (1, 4, 1, 1, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,12, 1, 6, 5), (1, 4, 1, 2, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,28, 1, 6, 5), (1, 4, 1, 6, 1))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1,  0, 1, 1, 5), (4, 1, 6, 1, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1, -4, 1, 1, 5), (4, 1, 6, 2, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,-20, 1, 1, 5), (4, 1, 6, 6, 3))
    
    a_equal(obj.convert_nbr_gij(ne, ngq,  0, 0, 1, 1, 5), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25, 0, 1, 1, 5), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  0,25, 1, 1, 5), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25,25, 1, 1, 5), (-1, -1, -1, -1, -1))

    # for paper
    ne = 3
    a_equal(obj.convert_nbr_gij(ne, ngq, 10, -11, 1, 2, 5), (3, 4, 1, 2, 3))



def test_convert_nbr_gij_6_6():
    '''
    cube_neighbor: convert_nbr_gij(): ne=6, np=4, panel=6
    '''
    ne = 6
    ngq = 4

    a_equal(obj.convert_nbr_gij(ne, ngq, gi=1, gj=1, ei=1, ej=3, panel=6), (1, 1, 1, 3, 6))
    a_equal(obj.convert_nbr_gij(ne, ngq, 4, 4, 6, 3, 6), (4, 4, 6, 3, 6))

    a_equal(obj.convert_nbr_gij(ne, ngq,  8, 1, 6, 3, 6), (1, 1, 3, 6, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 12, 1, 6, 3, 6), (1, 1, 3, 5, 2))
    a_equal(obj.convert_nbr_gij(ne, ngq, 28, 1, 6, 3, 6), (1, 1, 3, 1, 2))

    a_equal(obj.convert_nbr_gij(ne, ngq,   0, 1, 1, 3, 6), (4, 4, 4, 6, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq,  -4, 1, 1, 3, 6), (4, 4, 4, 5, 4))
    a_equal(obj.convert_nbr_gij(ne, ngq, -20, 1, 1, 3, 6), (4, 4, 4, 1, 4))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1, 8, 1, 6, 6), (4, 1, 6, 6, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,12, 1, 6, 6), (4, 1, 6, 5, 3))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,28, 1, 6, 6), (4, 1, 6, 1, 3))

    a_equal(obj.convert_nbr_gij(ne, ngq, 1,  0, 1, 1, 6), (1, 4, 1, 6, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1, -4, 1, 1, 6), (1, 4, 1, 5, 1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 1,-20, 1, 1, 6), (1, 4, 1, 1, 1))
    
    a_equal(obj.convert_nbr_gij(ne, ngq,  0, 0, 1, 1, 6), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25, 0, 1, 1, 6), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq,  0,25, 1, 1, 6), (-1, -1, -1, -1, -1))
    a_equal(obj.convert_nbr_gij(ne, ngq, 25,25, 1, 1, 6), (-1, -1, -1, -1, -1))
