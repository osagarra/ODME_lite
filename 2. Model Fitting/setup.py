# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:12:23 2013

@author: Oleguer Sagarra <osagarra@ub.edu> 
"""

from setuptools import setup, find_packages
readme = open('README.md').read()
setup(name='multi_edge_fitter',
      version='0.3',
      author='Oleguer Sagarra',
      author_email='osagarra@ub.edu',
      #license='',
      description="Module to solve the saddle point equaitons for maximum entropy multi-edge models",
      long_description=readme,
      packages=find_packages())
