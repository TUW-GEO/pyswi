#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for pyswi.

    This file was generated with PyScaffold 2.5.7, a tool that easily
    puts up a scaffold for your new Python project. Learn more under:
    http://pyscaffold.readthedocs.org/
"""

import sys
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
        cythonize(['pyswi/swi_img/swi_calc_routines.pyx'])


class NumpyBuildExt(_build_ext):

    def build_extensions(self):
        numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')

        for ext in self.extensions:
            if hasattr(ext, 'include_dirs') and not numpy_incl in ext.include_dirs:
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)

ext_modules = [Extension("pyswi.swi_img.swi_calc_routines",
                         ["pyswi/swi_img/swi_calc_routines.c"], include_dirs=[]), ]



def setup_package():
    cmdclass = {}
    cmdclass['cythonize'] = Cythonize
    cmdclass['build_ext'] = NumpyBuildExt
    needs_sphinx = {'build_sphinx', 'upload_docs'}.intersection(sys.argv)
    sphinx = ['sphinx'] if needs_sphinx else []
    setup(setup_requires=['six', 'pyscaffold>=2.5a0,<2.6a0'] + sphinx,
          tests_require=['test_swi'],
          cmdclass=cmdclass,
          ext_modules=ext_modules,
          use_pyscaffold=True)


if __name__ == "__main__":
    setup_package()
