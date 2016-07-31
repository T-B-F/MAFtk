"""Distutils based setup script for MAFtk.

Heavily inspired from the the Biopython setup.py script, all right reserved.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install

Or for more details about the options available from distutils, look at
the 'Installing Python Modules' distutils documentation, available from:

http://python.org/sigs/distutils-sig/doc/

MAFtk module is a set of function to be used with the MafIO parser
provided in biopython (supported in the branch MafIO of 
https://github.com/T-B-F/biopython)

This module contains read and write function to create a interval tree
based on the genomic position of the aligned blocks in MAF files.
A function is also provided to query the tree for a given species,
returning the alignments found from the query's start and stop.
The returned aligned are sliced according to the query's position.

Copyright 2016 by Tristan Bitard-Feildel. All rights reserved

Licensed under the MIT License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

https://opensource.org/licenses/MIT

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Dependencies
------------

MAFtk used the intervaltree library, copyright 2013-2015 Chaim-Leib Halbert,
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0
   
Biopython, licenced under the Biopython License Agreement 
(https://github.com/biopython/biopython/blob/master/LICENSE), should
also be installed

The branch MafIO that can be found at
https://github.com/T-B-F/biopython/tree/MafIO is required.


"""
from __future__ import print_function

import sys
import os
import shutil
import importlib

if "bdist_wheel" in sys.argv:
    try:
        import setuptools
        import wheel
    except ImportError:
        sys.exit("We need both setuptools AND wheel packages installed for bdist_wheel to work")
    # Import specific bits of setuptools ...
    from setuptools import setup
    from setuptools import Command
    from setuptools.command.install import install
    from setuptools.command.build_py import build_py
    from setuptools.command.build_ext import build_ext
    from setuptools import Extension
else:
    # Except for wheels, stick with standard library's distutils
    from distutils.core import setup
    from distutils.core import Command
    from distutils.command.install import install
    from distutils.command.build_py import build_py
    from distutils.command.build_ext import build_ext
    from distutils.extension import Extension


_CHECKED = None


def can_import(module_name):
    """can_import(module_name) -> module or None"""
    try:
        return importlib.import_module(module_name)
    except ImportError:
        return None


def is_intervaltree_installed():
    return bool(can_import("intervaltree"))

def is_biopython_installed():
    return bool(can_import("Bio"))

def is_mafio_installed():
    return bool(can_import("Bio.AlignIO.MafIO"))

_CHECKED = None

def is_automated():
    """Check for installation with easy_install or pip.
    """
    is_automated = False
    # easy_install: --dist-dir option passed
    try:
        dist_dir_i = sys.argv.index("--dist-dir")
    except ValueError:
        dist_dir_i = None
    if dist_dir_i is not None:
        dist_dir = sys.argv[dist_dir_i + 1]
        if "egg-dist-tmp" in dist_dir:
            is_automated = True
    # pip -- calls from python directly with "-c"
    if sys.argv in [["-c", "develop", "--no-deps"],
                    ["--no-deps", "-c", "develop"],
                    ["-c", "egg_info"]] \
                    or "pip-egg-info" in sys.argv \
                    or sys.argv[:3] == ["-c", "install", "--record"] \
                    or sys.argv[:4] == ['-c', 'install', '--single-version-externally-managed',
                                        '--record']:
        is_automated = True
    return is_automated

def check_dependencies_once():
    # Call check_dependencies, but cache the result for subsequent
    # calls.
    global _CHECKED
    if _CHECKED is None:
        _CHECKED = check_dependencies()
    return _CHECKED

def check_dependencies():
    
    if is_intervaltree_installed() and is_biopython_installed() and is_mafio_installed():
        return True
    if is_automated():
        return True
    
    if not is_intervaltree_installed():        
        print("""intervaltree is not installed

This package is required for MAFtk. Please install it before you install MAFtk. 

You can find intervaltree at https://github.com/chaimleib/intervaltree
or use pip install intervaltree
""")
    elif not is_biopython_installed():        
        print("""BioPython (Bio) is not installed

This package is required for MAFtk. Please install it before you install MAFtk. 

The biopython version required for the MafIO library can be find 
at https://github.com/T-B-F/biopython/tree/MafIO is required.
""")
    elif not is_mafio_installed():
        print("""The correct verion of BioPython (Bio) is not installed

This package is required for MAFtk. Please install it before you install MAFtk. 

The biopython version required for the MafIO library can be find 
at https://github.com/T-B-F/biopython/tree/MafIO is required.
""")
    sys.exit(1)


class install_maftk(install):
    """Override the standard install to check for dependencies.

    """
    user_options = install.user_options + [
        ('single-version-externally-managed', None,
            "used by system package builders to create 'flat' eggs"),
    ]
    boolean_options = install.boolean_options + [
        'single-version-externally-managed',
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.single_version_externally_managed = None

    def run(self):
        if check_dependencies():
            # Run the normal install.
            install.run(self)

class build_py_maftk(build_py):
    def run(self):
        if not check_dependencies_once():
            return
        build_py.run(self)


class build_ext_maftk(build_ext):
    def run(self):
        if not check_dependencies_once():
            return
        build_ext.run(self)


class test_maftk(Command):
    """Run the tests for the package.

    python setup.py build
    python setup.py install
    python setup.py test

    """
    description = "Automatically run the test suite for MAFtk."
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        
        os.chdir("tests")
        sys.path.insert(0, '')
        import run_tests
        run_tests.main([])

        # change back to the current directory
        os.chdir(this_dir)

# version is defined in lib/MAFtk/__init__.py
old_path = os.getcwd()

__version__ = "Undefined"
for line in open('lib/maftk/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

PACKAGES = ['maftk']

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = os.path.join('lib', *pkg.split('.'))


setup_args = {
    "name": 'MAFtk',
    "version": __version__,
    "author": 'Tristan Bitard-Feildel',
    "author_email": 'tristan@bitardfeildel.fr',
    "url": 'https://github.com/T-B-F/MAFtk',
    "description": 'Toolkit to manipulate genome alignments in MAF',
    "download_url": 'https://github.com/T-B-F/MAFtk/tarball/' + __version__,
    "cmdclass": {
        "install": install_maftk,
        "build_py": build_py_maftk,
        "build_ext": build_ext_maftk,
        "test": test_maftk,
        },
    "packages": PACKAGES,
    "package_dir": PACKAGE_DIR,
    "provides": ['maftk({0:s})'.format(__version__)],
    #"ext_modules": EXTENSIONS,
   }

try:
    setup(**setup_args)
finally:
    del sys.path[0]
    os.chdir(old_path)
