from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [#Extension("code_functions", ["code_functions.pyx"]),
                   Extension("powerline_code", ["powerline_code.pyx"]),
                   Extension("config", ["config.pyx"])]
)
