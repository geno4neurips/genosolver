from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np
import distutils.sysconfig
import platform

# eigenIncludeDir = "../eigen-eigen-5a0156e40feb"
eigenIncludeDir = "eigen"
# compilerFlags = "-std=c++11"

# cfg_vars = distutils.sysconfig.get_config_vars()
# for key, value in cfg_vars.items():
#     if type(value) == str:
#         cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


ext = Extension("genointerface",
                sources=["genointerface.pyx",
                         "pygenointerface.cpp",
                         "pygenonlp.cpp",
                         "lbfgsb.cpp",
                         "lineSearch.cpp",
                         "augmentedLagrangian.cpp"],
                language="c++",
                # extra_compile_args=[compilerFlags],
                include_dirs=['.',
                              np.get_include(),
                              # '..',
                              eigenIncludeDir]
                )
setup(name="genointerface",
      ext_modules=cythonize(ext))
