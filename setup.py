# -*- coding: utf-8 -*-
"""
    Setup file for pyswi.
    Use setup.cfg to configure your project.

    This file was generated with PyScaffold 3.2.3.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
import sys

from pkg_resources import VersionConflict, require
from setuptools import setup

from distutils.cmd import Command
from distutils.command.build_ext import build_ext as _build_ext
from distutils.extension import Extension
import pkg_resources


class Cythonize(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Make sure the compiled Cython files in the distribution are
        # up-to-date
        from Cython.Build import cythonize
        cythonize(['pyswi/swi_ts/swi_calc_routines.pyx'])


class NumpyBuildExt(_build_ext):

    def build_extensions(self):
        numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')

        for ext in self.extensions:
            if hasattr(ext, 'include_dirs') and not numpy_incl in ext.include_dirs:
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)


ext_modules = [Extension("pyswi.swi_ts.swi_calc_routines",
                         ["pyswi/swi_ts/swi_calc_routines.c"], include_dirs=[]), ]


try:
    require('setuptools>=38.3')
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)


if __name__ == "__main__":
    cmdclass = {}
    cmdclass['cythonize'] = Cythonize
    cmdclass['build_ext'] = NumpyBuildExt
    setup(cmdclass=cmdclass,
          ext_modules=ext_modules,
          use_pyscaffold=True)
