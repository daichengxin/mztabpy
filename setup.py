#!/usr/bin/env python
"""
The mztabpy library is a Python library that enables to handle of mztab files including conversion from tab-delimited format to
hdf5.
"""

from setuptools import setup, find_packages

version = '0.0.1'

def readme():
    with open('README.md') as f:
        return f.read()
setup(
    name='mztabpy',
    version=version,
    description='Python package to handle mztab files',
    author='Yasset Perez-Riverol',
    author_email='ypriverol@gmail.com',
    long_description=readme(),
    long_description_content_type='text/markdown',
    keywords='Proteomics, Label-free, quality control, MultiQC',
    url='https://github.com/bigbio/mztabpy/',
    download_url='https://github.com/bigbio/mztabpy/',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'tables',
        'pandas',
        'Click'
    ],
    scripts=['mztabpy_click.py'],
    entry_points={
        'console_scripts': [
          'mztabpy = mztabpy.mztabpy_cli:main'
        ]},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)