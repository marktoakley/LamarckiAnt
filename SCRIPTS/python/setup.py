import os
import numpy as np
from distutils.core import setup
from distutils.extension import Extension

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

setup(
    name="ptmc",
    version='0.0.1', 
    description="tools for parallel tempering monte carlo with GMIN",
    #cmdclass = {'build_ext': build_ext},
    packages=["ptmc",
              "ptmc.src",
             ],
    ext_modules= 
        [
            Extension("ptmc.src.exch_tools", ["ptmc/src/exch_tools.c", "ptmc/src/_exch_tools.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O2',],
                        ),
        ]
      )
