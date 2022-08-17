from setuptools import setup
from Cython.Build import cythonize
import numpy
import sys

sys.setrecursionlimit(20000)

setup(
    name='Skyrme Polynomial',
    ext_modules=cythonize("skyrme_poly_2.pyx"),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
    compiler_directives={'language_level' : "3",'extra_compile_args' : "-O3"},
)