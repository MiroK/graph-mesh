#!/usr/bin/env python

from setuptools import find_packages, setup

setup(name = 'graph_mesh',
      version = '0.1',
      description = 'FEniCS meshes from graphs',
      author = 'Miroslav Kuchta',
      author_email = 'miroslav.kuchta@gmail.com',
      url = 'https://github.com/mirok/graph-mesh.git',
      packages=find_packages(),
      include_package_data=True)
