# QMCkl: Quantum Monte Carlo Kernel Library

![Build Status](https://github.com/TREX-CoE/qmckl/workflows/test-build/badge.svg?branch=main)

The domain of quantum chemistry needs a library in which the main
kernels of Quantum Monte Carlo (QMC) methods are implemented. In the
library proposed in this project, we expose the main algorithms in a
simple language and provide a standard API and tests to enable the
development of high-performance QMCkl implementations taking
advantage of modern hardware.

See the [source code](https://github.com/TREX-CoE/qmckl/tree/main/src)
to read the documentation.


To clone the repository, use:
```
git clone https://github.com/TREX-CoE/qmckl.git
```

# Installation

The simplest way to obtain the source files of QMCkl is to download a source
distribution. This particular repository is for maintainers, who write the kernels
in org-mode files and produce the source code and the documentation from these files.

## For maintainers

```
./autogen.sh
QMCKL_DEVEL=1 ./configure --prefix=$PWD/_install --enable-silent-rules --enable-maintainer-mode

make
make check
```

## For users

Obtain a source distribution and run

```
./configure 
make
make check
sudo make install
sudo make installcheck
```

------------------------------

![European flag](https://trex-coe.eu/sites/default/files/inline-images/euflag.jpg)
[TREX: Targeting Real Chemical Accuracy at the Exascale](https://trex-coe.eu) project has received funding from the European Union’s Horizon 2020 - Research and Innovation program - under grant agreement no. 952165. The content of this document does not represent the opinion of the European Union, and the European Union is not responsible for any use that might be made of such content.

