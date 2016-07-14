from distutils.core import setup
from Cython.Build import cythonize

import os
currentDir = os.path.dirname(os.path.realpath(__file__))

setup(
    ext_modules = cythonize(os.path.join(currentDir, "helloworld.pyx"))
)