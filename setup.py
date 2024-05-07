#!/usr/bin/env python

from setuptools import setup
exec(open('idattools/__init__.py').read())


setup(
    name='idat-tools',
    scripts=['bin/idat-tools'],
    packages=["idattools"],
    version=__version__,
    author=__author__,
    url=__homepage__,
    description='Toolkit to read, modify and export idat files',
    long_description=open("README.md", 'r').read().strip(),
    setup_requires=['setuptools'],# bit odd, this can only be loaded if it is there
    install_requires=[_.strip() for _ in open("requirements.txt", "r").readlines() if _[0] != "#"],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)

