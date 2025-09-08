#!/usr/bin/env python

import os
import re
from setuptools import setup
import subprocess

"""
23.04.20 add alphashape
"""


setup(name='CellModeller',
    install_requires=['numpy<2', 'scipy', 'pyopengl', 'mako', 'pyqt5', 'pyopencl', 'reportlab', 'matplotlib',
                      'vtk',
                      'pyabc',
                      'opencv-python<4.10',
                      'scikit-image',
                      'alphashape'],
    setup_requires=['numpy<2', 'scipy', 'pyopengl', 'mako', 'pyqt5', 'pyopencl', 'reportlab', 'matplotlib',
                    'vtk',
                    'pyabc',
                    'opencv-python<4.10',
                    'scikit-image',
                    'alphashape'],
    packages=['CellModeller',
                'CellModeller.Biophysics',
                'CellModeller.Biophysics.BacterialModels',
                'CellModeller.Biophysics.GeneralModels',
                'CellModeller.Integration',
                'CellModeller.Regulation',
                'CellModeller.Signalling',
                'CellModeller.GUI'],
     package_data={'':['*.cl','*.ui']},
     python_requires='>=3',
     version="0.0.1"
)
