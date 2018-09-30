#!/usr/bin/env python

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='',top_path=None):
    config = Configuration('',parent_package,top_path)
    config.add_extension('c_routines',
                         sources = ['src/c/c_routines.pyf','src/c/c_routines.c'])
    return config

if __name__ == "__main__":
    setup(**configuration(top_path='').todict())
