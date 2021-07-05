#!/bin/bash
# Script to build the final qmckl.h file

# All the produced header files are concatenated in the =qmckl.h=
# file, located in the include directory. The =*_private.h= files
# are excluded.

set -e


# Check required environment variables
# ------------------------------------

if [[ -z ${srcdir} ]] ; then
   echo "Error: srcdir environment variable is not defined"
   exit 1
fi

if [[ -z ${qmckl_h} ]] ; then
   echo "Error: qmckl_h environment variable is not defined"
   exit 1
fi


# Generate C header file
# ----------------------

HEADERS=""
for i in $(cat ${srcdir}/org/table_of_contents)
do
    HEADERS+="${i%.org}_type.h "
done

for i in $(cat ${srcdir}/org/table_of_contents)
do
    HEADERS+="${i%.org}_func.h "
done

OUTPUT=${qmckl_h}

cat << EOF > ${OUTPUT}
/*
 *    ------------------------------------------
 *     QMCkl - Quantum Monte Carlo kernel library
 *     ------------------------------------------
 *
 *     Documentation : https://trex-coe.github.io/qmckl
 *     Issues        : https://github.com/trex-coe/qmckl/issues
 *
 *     BSD 3-Clause License
 *
 *     Copyright (c) 2020, TREX Center of Excellence
 *     All rights reserved.
 *
 *     Redistribution and use in source and binary forms, with or without
 *     modification, are permitted provided that the following conditions are met:
 *
 *     1. Redistributions of source code must retain the above copyright notice, this
 *        list of conditions and the following disclaimer.
 *
 *     2. Redistributions in binary form must reproduce the above copyright notice,
 *        this list of conditions and the following disclaimer in the documentation
 *        and/or other materials provided with the distribution.
 *
 *     3. Neither the name of the copyright holder nor the names of its
 *        contributors may be used to endorse or promote products derived from
 *        this software without specific prior written permission.
 *
 *     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 *
 */

#ifndef __QMCKL_H__
#define __QMCKL_H__

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
EOF

for i in ${HEADERS}
do
    header=${srcdir}/src/$i
    if [[ -f $header ]] ; then
        echo "/* $header */" >> ${OUTPUT}
        cat $header >> ${OUTPUT}
    fi
done

cat << EOF >> ${OUTPUT}
#endif
EOF


