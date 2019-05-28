# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np
# import distutils.sysconfig
# import platform
import ssl
import wget
import os
import zipfile
import sys

if sys.platform.startswith("linux"):
    pass
    # os.environ["CC"] = "clang++"
    # os.environ["CXX"] = "clang++"
elif sys.platform.startswith("win"):
    pass
elif sys.platform.startswith("darwin"):
    pass
else:
    # TODO do not set any flags,
    # recommend to setup their own compilation script,
    # if this one fails
    raise Exception("OS not supported!")

# Notes
# clang works only with "-std=c++11"
# what if we are on windows ???

dir_eigen_include = os.path.join("genosolver", "eigen")
url = "https://github.com/eigenteam/eigen-git-mirror/archive/3.3.7.zip"
file_eigen_zip = "eigen-git-mirror-3.3.7.zip"

# download eigen if no path is provided
if not os.path.isdir(dir_eigen_include):
    if not os.path.isfile(file_eigen_zip):
        print("Downloading eigen requirements...")
        ssl._create_default_https_context = ssl._create_unverified_context
        file_eigen_zip = wget.download(url, bar=None)
        print("\n done.")
    print("Extracting files")
    with zipfile.ZipFile(file_eigen_zip, "r") as zip_file:
        zip_file.extractall("genosolver")
    os.rename(os.path.join("genosolver", file_eigen_zip[:-4]),
              dir_eigen_include)
    os.remove(file_eigen_zip)


ext = Extension(name="genosolver.genointerface",
                sources=[
                    os.path.join("genosolver", "genointerface.pyx"),
                    os.path.join("genosolver", "pygenointerface.cpp"),
                    os.path.join("genosolver", "pygenonlp.cpp"),
                    os.path.join("genosolver", "lbfgsb.cpp"),
                    os.path.join("genosolver", "lineSearch.cpp"),
                    os.path.join("genosolver", "augmentedLagrangian.cpp")],
                language="c++",
                extra_compile_args=["-std=c++11"],
                include_dirs=["genosolver",
                              np.get_include(),
                              dir_eigen_include]
                )
setup(name="genosolver",
      version="1.0",
      description="GENO (GENeric Optimization) solver",
      url="http://geno4neurips.pythonanywhere.com",
      author="GENO team",
      license="",
      zip_safe=False,
      packages=find_packages(),
      ext_modules=cythonize(ext),
      python_requires=">3.7")

os.remove(os.path.join("genosolver", "genointerface.cpp"))

# # from distutils.core import setup, Extension
# from setuptools import setup, Extension, find_packages
# from Cython.Build import cythonize
# import numpy as np
# import os

# # import distutils.sysconfig
# # import platform

# # dir_eigen_inc = os.path.join("genosolver", "eigen")


# # os.chdir("genosolver")
# ext = Extension("genosolver.genointerface",
#                 sources=[
#                     "genosolver/pygenointerface.cpp",
#                     "genosolver/pygenonlp.cpp",
#                     "genosolver/lbfgsb.cpp",
#                     "genosolver/lineSearch.cpp",
#                     "genosolver/augmentedLagrangian.cpp",
#                     "genosolver/genointerface.pyx"],
#                 language="c++",
#                 # extra_compile_args=[compilerFlags],
#                 include_dirs=[np.get_include(),
#                               ".",
#                               "genosolver/genosolver"
#                               "genosolver/*"
#                               "genosolver",
#                               "genosolver/eigen"]
#                 )
# # os.chdir("..")
# setup(name="genosolver",
#       packages=find_packages(),
#       ext_modules=cythonize(ext))


# os.remove("genosolver/genointerface.cpp")
