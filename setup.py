#!/usr/bin/env python
from distutils.core import setup

with open('README.md') as file:
	long_description = file.read()

setup(name = 'dm2linear',
      version = '0.1',
      description = 'Expand density to linear combination of GTO from density matrix',
      long_description = long_description,
      author = 'Li Li',
      author_email = 'aiki.nogard@gmail.com',
      install_requires = ['numpy'],
      packages = ['dm2linear']
      ) 