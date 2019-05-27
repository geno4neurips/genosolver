# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np
import os

# import distutils.sysconfig
# import platform

dir_eigen_inc = os.path.join("eigen")

ext = Extension("genointerface",
                sources=[
                    "genointerface.pyx",
                    "pygenointerface.cpp",
                    "pygenonlp.cpp",
                    "lbfgsb.cpp",
                    "lineSearch.cpp",
                    "augmentedLagrangian.cpp"],
                language="c++",
                # extra_compile_args=[compilerFlags],
                include_dirs=[np.get_include(),
                              ".",
                              # "genosolver/"
                              # "genosolver",
                              dir_eigen_inc]
                )
setup(name="genosolver",
      # packages=find_packages(),
      ext_modules=cythonize(ext))

os.remove("genointerface.cpp")
