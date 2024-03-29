'''

abstract : call the subroutines and functions from a f90 library

history :
  2017-06-29  ki-hwan kim  split from f90 wrapper files

'''

from __future__ import print_function, division
from ctypes import POINTER, c_int, c_bool, c_double, c_float
from ctypes import CDLL, RTLD_GLOBAL
import logging
import numpy as np
import numpy.ctypeslib as npct




def fmod2py(so_dpath, libname, modname, func_args, fdtype='f8'):
    '''
    wrapper of a f90 module to call subroutines and functions
    '''

    f90_lib = npct.load_library(libname, so_dpath)
    #f90_lib = ctypes.CDLL(so_dpath+'/'+libname+'.so', mod=RTLD_GLOBAL)

    ptrs = {'i': POINTER(c_int),
            'bool': POINTER(c_bool),
            'i1d': npct.ndpointer(ndim=1, dtype='i4'),
            'i2d': npct.ndpointer(ndim=2, dtype='i4', flags='F'),
            'i3d': npct.ndpointer(ndim=3, dtype='i4', flags='F'),
            'i4d': npct.ndpointer(ndim=4, dtype='i4', flags='F'),
            'i5d': npct.ndpointer(ndim=5, dtype='i4', flags='F'),
            'f1d': npct.ndpointer(ndim=1, dtype=fdtype, flags='F'),
            'f2d': npct.ndpointer(ndim=2, dtype=fdtype, flags='F'),
            'f3d': npct.ndpointer(ndim=3, dtype=fdtype, flags='F'),
            'f4d': npct.ndpointer(ndim=4, dtype=fdtype, flags='F'),
            'f5d': npct.ndpointer(ndim=5, dtype=fdtype, flags='F'),
            'f6d': npct.ndpointer(ndim=6, dtype=fdtype, flags='F')}

    if fdtype in [np.float64, 'f8']:
        ptrs['f'] = POINTER(c_double)
    elif fdtype in [np.float32, 'f4']:
        ptrs['f'] = POINTER(c_float)
    else:
        raise ValueError('The fdtype is invalid: {}'.format(fdtype))

    func_dict = dict()
    for funcname, (argtypes, restype) in func_args.items():
        func = getattr(f90_lib, '__{}_MOD_{}'.format(modname, funcname))

        if argtypes == None:
            func.argtypes = None
        else:
            func.argtypes = [ptrs[argtype] for argtype in argtypes]

        func.restype = {
                None: None,
                'i': c_int,
                'f': c_double if fdtype in [np.float64, 'f8'] else c_float,
                'bool': c_bool}[restype]

        func_dict[funcname] = func

    return func_dict
