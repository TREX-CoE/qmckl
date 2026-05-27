#!/usr/bin/env python3
"""
setup.py file for qmckl package (pure Python, ctypes-based).
"""

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="qmckl",
    version="1.1.0",
    author="TREX-CoE",
    author_email="posenitskiy@irsamc.ups-tlse.fr",
    description="Python API of the QMCkl library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    url="https://github.com/TREX-CoE/qmckl",
    license="BSD",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS",
    ],
    python_requires=">=3.6",
    install_requires=["numpy"],
)
