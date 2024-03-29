'''

abstract : wrapper to cube_neighbor.f90

history :
  2017-04-15  ki-hwan kim  start
  2017-05-09  ki-hwan kim  remove convert_neighbor_ij(), find_nbr_elem_ij()
                           add convert_nbr_elem_ij()
  2017-05-15  ki-hwan kim  convert_nbr_elem_ij() -> convert_nbr_eij()
                           add convert_nbr_gij()

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join
from ctypes import c_int, byref
import logging
import numpy as np

from f90wrap import fmod2py




class CubeNeighbor(object):
    ''' 
    Wrapper the cube_neighbor.f90 library
    '''

    def __init__(self):
        self.ngq = 4
        
        #
        # Load the library using the numpy ctypeslib
        #
        so_dpath = join(dirname(abspath(__file__)), 'f90')
        libname = 'libshared'
        modname = 'cube_neighbor'
        func_args = { \
            'quotient': [('i','i'), 'i'],
            'init_nbr_panels': [None, None],
            'convert_rotated_ij': [('i','i','i','i','i1d'), None],
            'convert_nbr_eij': [('i','i','i','i','i1d'), None],
            'convert_nbr_gij': [('i','i','i','i','i','i','i','i1d'), None]}

        self.f90_funcs = fmod2py(so_dpath, libname, modname, func_args)

        self.f90_funcs['init_nbr_panels']()


    def quotient(self, n, i):
        n_p = byref(c_int(n))
        i_p = byref(c_int(i))

        q = self.f90_funcs['quotient'](n_p, i_p)

        return q


    def convert_rotated_ij(self, n, i, j, rot):
        n_p = byref(c_int(n))
        i_p = byref(c_int(i))
        j_p = byref(c_int(j))
        rot_p = byref(c_int(rot))
        ret = np.zeros(2, 'i4')  # rotated (i,j)
        self.f90_funcs['convert_rotated_ij'](n_p, i_p, j_p, rot_p, ret)

        return ret


    def convert_nbr_eij(self, ne, ei, ej, panel):
        ne_p = byref(c_int(ne))
        ei_p = byref(c_int(ei))
        ej_p = byref(c_int(ej))
        p_p  = byref(c_int(panel))
        ret = np.zeros(4, 'i4')  # rotated elem (ei, ej, panel, rot)
        self.f90_funcs['convert_nbr_eij'](ne_p, ei_p, ej_p, p_p, ret)

        return ret


    def convert_nbr_gij(self, ne, ngq, gi, gj, ei, ej, panel):
        ne_p = byref(c_int(ne))
        ngq_p = byref(c_int(ngq))
        gi_p = byref(c_int(gi))
        gj_p = byref(c_int(gj))
        ei_p = byref(c_int(ei))
        ej_p = byref(c_int(ej))
        p_p  = byref(c_int(panel))
        ret = np.zeros(5, 'i4')  # rotated point (gi, gj, ei, ej, panel)
        self.f90_funcs['convert_nbr_gij'](ne_p, ngq_p, gi_p, gj_p, ei_p, ej_p, p_p, ret)

        return ret
