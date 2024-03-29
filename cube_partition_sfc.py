'''

abstract : wrapper to cube_partition_sfc.f90

history :
  2017-04-06  ki-hwan kim  start
  2017-06-16  ki-hwan kim  add make_elem_coord()

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join
from ctypes import c_int, byref
import logging
import numpy as np

from f90wrap import fmod2py




class CubePartitionSFC(object):
    ''' 
    Wrapper the cube_partition_sfc.f90 library
    '''

    def __init__(self, ne, nproc):
        self.ne = ne
        self.nproc = nproc

        #
        # Load the library using the numpy ctypeslib
        #
        so_dpath = join(dirname(abspath(__file__)), 'f90')
        libname = 'libshared'
        modname = 'cube_partition_sfc'
        func_args = { \
            'rot': [('i','i','i2d','i2d'), None],
            'inv_x': [('i','i2d','i2d'), None],
            'inv_y': [('i','i2d','i2d'), None],
            'make_sfcs': [('i3d','i3d','i3d'), None],
            'find_size_factors': [('i',), 'i'],
            'find_factors': [('i','i','i','i1d'), None],
            'make_panel_sfc': [('i','i','i2d'), None],
            'make_global_sfc': [('i','i','i3d'), None],
            'make_cube_rank': [('i','i','i1d','i3d','i3d'), None],
            'make_elem_coord': [('i','i','i','i3d','i3d','i2d'), None]}

        self.f90_funcs = fmod2py(so_dpath, libname, modname, func_args)


    def rot(self, num_rot, arr):
        n_p = byref(c_int(arr.shape[0]))
        num_rot_p = byref(c_int(num_rot))
        ret = np.zeros(arr.shape, 'i4', order='F')
        self.f90_funcs['rot'](n_p, num_rot_p, arr, ret)

        return ret


    def inv_x(self, arr):
        n_p = byref(c_int(arr.shape[0]))
        ret = np.zeros(arr.shape, 'i4', order='F')
        self.f90_funcs['inv_x'](n_p, arr, ret)

        return ret


    def inv_y(self, arr):
        n_p = byref(c_int(arr.shape[0]))
        ret = np.zeros(arr.shape, 'i4', order='F')
        self.f90_funcs['inv_y'](n_p, arr, ret)

        return ret


    def make_sfcs(self):
        hilbert = np.zeros((2,2,4), 'i4', order='F')
        peano = np.zeros((3,3,4), 'i4', order='F')
        cinco = np.zeros((5,5,4), 'i4', order='F')
        self.f90_funcs['make_sfcs'](hilbert, peano, cinco)

        return hilbert, peano, cinco


    def find_size_factors(self):
        ne = self.ne
        ne_p = byref(c_int(ne))
        return self.f90_funcs['find_size_factors'](ne_p)


    def find_factors(self):
        ne = self.ne
        nproc = self.nproc
        size_factors = self.find_size_factors()
        ne_p = byref(c_int(ne))
        nproc_p = byref(c_int(nproc))
        size_p = byref(c_int(size_factors))
        ret = np.zeros(size_factors, 'i4', order='F')
        self.f90_funcs['find_factors'](ne_p, nproc_p, size_p, ret)

        return ret


    def make_panel_sfc(self):
        ne = self.ne
        nproc = self.nproc
        ne_p = byref(c_int(ne))
        nproc_p = byref(c_int(nproc))
        ret = np.zeros((ne,ne), 'i4', order='F')
        self.f90_funcs['make_panel_sfc'](ne_p, nproc_p, ret)

        return ret


    def make_global_sfc(self):
        ne = self.ne
        nproc = self.nproc
        ne_p = byref(c_int(ne))
        nproc_p = byref(c_int(nproc))
        ret = np.zeros((ne,ne,6), 'i4', order='F')
        self.f90_funcs['make_global_sfc'](ne_p, nproc_p, ret)

        return ret


    def make_cube_rank(self):
        ne = self.ne
        nproc = self.nproc
        ne_p = byref(c_int(ne))
        nproc_p = byref(c_int(nproc))
        nelems = np.zeros(nproc, 'i4', order='F')
        cube_rank = np.zeros((ne,ne,6), 'i4', order='F')
        cube_lid = np.zeros((ne,ne,6), 'i4', order='F')
        self.f90_funcs['make_cube_rank'](ne_p, nproc_p, nelems, cube_rank, cube_lid)

        return nelems, cube_rank, cube_lid


    def make_elem_coord(self, iproc, nelem, cube_rank, cube_lid):
        '''
        return (ei,ej,panel)x(nelem)
        '''
        ne = self.ne

        ne_p = byref(c_int(ne))
        iproc_p = byref(c_int(iproc))
        nelem_p = byref(c_int(nelem))
        elem_coord = np.zeros((3,nelem), 'i4', order='F')
        self.f90_funcs['make_elem_coord'](ne_p, iproc_p, nelem_p, cube_rank, cube_lid, elem_coord)

        return elem_coord
