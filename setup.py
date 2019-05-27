# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np
import os

# import distutils.sysconfig
# import platform

# dir_eigen_inc = os.path.join("genosolver", "eigen")


# os.chdir("genosolver")
ext = Extension("genosolver.genointerface",
                sources=[
                    "genosolver/pygenointerface.cpp",
                    "genosolver/pygenonlp.cpp",
                    "genosolver/lbfgsb.cpp",
                    "genosolver/lineSearch.cpp",
                    "genosolver/augmentedLagrangian.cpp",
                    "genosolver/genointerface.pyx"],
                language="c++",
                # extra_compile_args=[compilerFlags],
                include_dirs=[np.get_include(),
                              ".",
                              "genosolver/genosolver"
                              "genosolver/*"
                              "genosolver",
                              "genosolver/eigen"]
                )
# os.chdir("..")
setup(name="genosolver",
      packages=find_packages(),
      ext_modules=cythonize(ext))


os.remove("genosolver/genointerface.cpp")
