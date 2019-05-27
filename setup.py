# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np
import os

# import distutils.sysconfig
# import platform

# eigenIncludeDir = "../eigen-eigen-5a0156e40feb"
dir_eigen_inc = os.path.join("genosolver", "eigen")
# compilerFlags = "-std=c++11"

# cfg_vars = distutils.sysconfig.get_config_vars()
# for key, value in cfg_vars.items():
#     if type(value) == str:
#         cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


ext = Extension("genointerface",
                sources=[
                    os.path.join("genosolver", "genointerface.pyx"),
                    os.path.join("genosolver", "pygenointerface.cpp"),
                    os.path.join("genosolver", "pygenonlp.cpp"),
                    os.path.join("genosolver", "lbfgsb.cpp"),
                    os.path.join("genosolver", "lineSearch.cpp"),
                    os.path.join("genosolver", "augmentedLagrangian.cpp")],
                language="c++",
                # extra_compile_args=[compilerFlags],
                include_dirs=[np.get_include(),
                              "genosolver",
                              dir_eigen_inc]
                )
setup(name="genosolver",
      packages=find_packages(),
      ext_modules=cythonize(ext))
