
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension(
              "cy_modules.cython_helper"     ,  # location of the resulting .so
              ["cy_modules/cython_helper.pyx"], # source file(s)
              libraries=["m"],
              #extra_compile_args = ["-O3","-ffast-math"] # dangerous
              extra_compile_args = ["-O3"],
             ),
    ]

setup(name        = 'cython_module',
      packages    = find_packages(),
      cmdclass    = {'build_ext': build_ext},
      ext_modules = ext_modules,
     )

# run with:
# python setup.py build_ext --inplace
# in case of strange errors:
# export CFLAGS="-I /usr/local/lib/python2.7/dist-packages/numpy/core/include $CGLAGS"
# use as:
# import cython_helper as bla
