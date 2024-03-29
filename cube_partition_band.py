'''

abstract : wrapper to cube_partition_band.f90

history :
  2018-03-06  ki-hwan kim  start

'''

from __future__ import print_function, division
from os.path import dirname, abspath, join
from ctypes import c_int, byref
import logging
import numpy as np

from f90wrap import fmod2py
from cube_neighbor import CubeNeighbor




class CubePartitionBand(object):
    ''' 
    Wrapper the cube_partition.f90 library
    '''

    def __init__(self, ne, nproc):
        self.ne = ne
        self.nproc = nproc

        #
        # Load the library using the numpy ctypeslib
        #
        so_dpath = join(dirname(abspath(__file__)), 'f90')
        libname = 'libshared'
        modname = 'cube_partition_band'
        func_args = { \
            'calc_perimeter_ratio': [('i','i','i','i','i','i1d','i1d','i2d'), 'f'],
            'find_optimal_band': [('i','i','i','i','i','i1d','i2d','i1d'), None],
            'band_partition': [('i','i','i1d','i3d'), None],
            'make_cube_rank': [('i','i','i1d','i3d','i3d'), None],
            'global_perimeter_ratio': [('i','i','i3d','i2d'), 'f'],
            'global_communication_ratio': [('i','i','i','i3d','i2d'), 'f'],
            'make_cube_color': [('i','i','i3d','i3d'), None]}

        self.f90_funcs = fmod2py(so_dpath, libname, modname, func_args)

        CubeNeighbor() # to initialize the cube_neighbor module


    def calc_perimeter_ratio(self, start_rank, end_rank, nelems, i12, box):
        ne = self.ne
        nproc = self.nproc
        nx, ny = box.shape

        to_i = lambda x: byref(c_int(x))
        to_f = lambda x: byref(c_double(x))

        perimeter_ratio = self.f90_funcs['calc_perimeter_ratio'](
                to_i(nx), to_i(ny), to_i(nproc), to_i(start_rank), to_i(end_rank),
                nelems, i12, box)

        return perimeter_ratio


    def find_optimal_band(self, start_rank, start_i, nelems, box):
        nx, ny = box.shape
        nproc = nelems.size

        to_i = lambda x: byref(c_int(x))

        ret = np.zeros(2, 'i4')
        self.f90_funcs['find_optimal_band'](
                to_i(nx), to_i(ny), to_i(nproc), to_i(start_rank), to_i(start_i),
                nelems, box, ret)

        return ret[0], ret[1]  # (rank, i2)


    def band_partition(self, nelems):
        ne = self.ne
        nproc = self.nproc

        to_i = lambda x: byref(c_int(x))

        cube_rank = np.zeros((ne,ne,6), 'i4', order='F')
        self.f90_funcs['band_partition'](
                to_i(ne), to_i(nproc), nelems, cube_rank)

        return cube_rank


    def make_cube_rank(self):
        ne = self.ne
        nproc = self.nproc

        to_i = lambda x: byref(c_int(x))

        nelems = np.zeros(nproc, 'i4')
        cube_rank = np.zeros((ne,ne,6), 'i4', order='F')
        cube_lid = np.zeros((ne,ne,6), 'i4', order='F')
        self.f90_funcs['make_cube_rank'](
                to_i(ne), to_i(nproc), nelems, cube_rank, cube_lid)

        return nelems, cube_rank, cube_lid


    def global_perimeter_ratio(self, cube_rank):
        ne = self.ne
        nproc = self.nproc

        to_i = lambda x: byref(c_int(x))

        num_nbrs = np.zeros((2,nproc), 'i4', order='F')
        perimeter_ratio = self.f90_funcs['global_perimeter_ratio'](
                to_i(ne), to_i(nproc), cube_rank, num_nbrs)

        return perimeter_ratio, num_nbrs



    def global_communication_ratio(self, ngq, cube_rank):
        ne = self.ne
        nproc = self.nproc

        to_i = lambda x: byref(c_int(x))

        num_pts = np.zeros((2,nproc), 'i4', order='F')
        comm_ratio = self.f90_funcs['global_communication_ratio'](
                to_i(ne), to_i(ngq), to_i(nproc), cube_rank, num_pts)

        return comm_ratio, num_pts



    def make_cube_color(self, cube_rank):
        ne = self.ne
        nproc = self.nproc

        to_i = lambda x: byref(c_int(x))

        cube_color = np.zeros((ne,ne,6), 'i4', order='F')
        perimeter_ratio = self.f90_funcs['make_cube_color'](
                to_i(ne), to_i(nproc), cube_rank, cube_color)

        return cube_color
