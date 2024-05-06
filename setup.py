#!/usr/bin/env python

from setuptools import setup


with open("README.md", 'r') as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as fh:
    requirements = [_.strip() for _ in fh.readlines() if _[0] != "#"]


setup(
    name='idat-tools',
    scripts=['bin/idat-tools'],
    packages=["idattools"],
    version='0.0',
    description='Toolkit to read, modify and export idat files',
    long_description=long_description,
    author='Youri Hoogstrate',
    url="https://github.com/yhoogstrate/idat-tools",
    install_requires=['wheel']
)
