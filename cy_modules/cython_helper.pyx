import cython
import numpy as N
cimport numpy as N
cimport cython

from libc.math cimport exp

#-------------------
# type definitions -
#-------------------

DTYPEf = N.float64
DTYPEl = N.int_

ctypedef N.float64_t DTYPEf_t
ctypedef N.int_t     DTYPEl_t


cpdef img_from_vis(
                    N.ndarray[DTYPEf_t, ndim=1] u,
                    N.ndarray[DTYPEf_t, ndim=1] v,
                    N.ndarray[DTYPEf_t, ndim=1] weight,
                    N.ndarray[DTYPEf_t, ndim=1] vis,
                    long                        imsize,
                    DTYPEf_t                    pixsize
                  ):
    """
      u  (nvis) : in meters
      v  (nvis) : in meters
      weight  (nvis) or (2 x nvis) : length is either twice that of u, v, with real_1, image_1, real_2, imag_2, real_3, ... or same as that of u, v
      vis  (nvis) or (2 x nvis) : length is either twice that of u, v, with real_1, image_1, real_2, imag_2, real_3, ... or same as that of u, v
      imsize
      pixsize

    """

    cdef DTYPEf_t arcsec2rad

    arcsec2rad = N.pi / 180./ 3600.
    xxx = N.linspace(-imsize/2 * pixsize * arcsec2rad,
                        imsize/2 * pixsize * arcsec2rad,
                        imsize)
    X, Y = N.meshgrid(xxx, xxx)

    # some CASA pixel alignment
    X -= pixsize / 2. * arcsec2rad
    Y += pixsize / 2. * arcsec2rad

    cdef N.ndarray[N.complex64_t, ndim=2] zz = N.zeros(X.shape).astype(N.complex64)
    cdef N.ndarray[DTYPEf_t, ndim=1] real = N.zeros((len(u)), dtype=DTYPEf)
    cdef N.ndarray[DTYPEf_t, ndim=1] imag = N.zeros((len(u)), dtype=DTYPEf)

    real = vis[0::2]
    imag = vis[1::2]

    if len(weight) == len(u):
        for ii in range(len(u)):
            zz += weight[ii] * (real[ii] + 1j * imag[ii]) * N.exp(2j * N.pi * (u[ii] * X + v[ii] * Y))
    else:
        # weight in format of real_1, image_1, real_2, imag_2, real_3 ....
        for ii in range(len(u)):
            zz += weight[::2][ii] * (real[ii] + 1j * imag[ii]) * N.exp(2j * N.pi * (u[ii] * X + v[ii] * Y))

    return N.array(zz, dtype=N.complex64)