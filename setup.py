
import os

from setuptools import setup

setup (
    name='zmattools',
    version='1.0.0',
    description='Tools to generate and manipulate z-matrices',
    long_description=open('README.md').read(),
    url='https://github.com/sgenheden/zmattools',
    author='Samuel Genheden',
    author_email='samuel.genheden@gmail.com',
    license='MIT',

    packages=['zmat',],
    entry_points={'console_scripts': ['zmat = zmat.tools:main_zmat', ]},
    install_requires=['numpy','parmed','networkx'],
)
