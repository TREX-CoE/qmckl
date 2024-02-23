#!/bin/sh
# Script to build the final src/qmckl_f.F90 file

set -e

# All the produced header files are concatenated in the =src/qmckl_f.F90=
# file, located in the share/qmckl/fortran directory.


# Check required environment variables
# ------------------------------------

if [ -z ${srcdir} ] ; then
   echo "Error: srcdir environment variable is not defined"
   exit 1
fi

if [ -z ${top_builddir} ] ; then
   echo "Error: top_builddir environment variable is not defined"
   exit 1
fi

if [ -z ${qmckl_f} ] ; then
   echo "Error: qmckl_f environment variable is not defined"
   exit 1
fi



# Generate Fortran interface file
# -------------------------------

HEADERS_TYPE="src/qmckl_*_fh_type.F90"
HEADERS="src/qmckl_*_fh_func.F90"

cat << EOF > ${qmckl_f}
!
!    ------------------------------------------
!     QMCkl - Quantum Monte Carlo kernel library
!     ------------------------------------------
!    
!     Documentation : https://trex-coe.github.io/qmckl
!     Issues        : https://github.com/trex-coe/qmckl/issues
!    
!     BSD 3-Clause License
!     
!     Copyright (c) 2020, TREX Center of Excellence
!     All rights reserved.
!     
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are met:
!     
!     1. Redistributions of source code must retain the above copyright notice, this
!        list of conditions and the following disclaimer.
!     
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     
!     3. Neither the name of the copyright holder nor the names of its
!        contributors may be used to endorse or promote products derived from
!        this software without specific prior written permission.
!     
!     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     
!     
!    
!    
!
module qmckl_constants
  use, intrinsic :: iso_c_binding
EOF

for i in ${HEADERS_TYPE}
do
    cat $i >> ${qmckl_f}
done

cat << EOF >> ${qmckl_f}
end module qmckl_constants

module qmckl
  use qmckl_constants
EOF

for i in ${HEADERS}
do
    cat $i >> ${qmckl_f}
done

cat << EOF >> ${qmckl_f}
end module qmckl
EOF
