# QMCkl - Quantum Monte Carlo kernel library
#
# BSD 3-Clause License
#
# Copyright (c) 2020, TREX Center of Excellence
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

dist_src_DATA = $(ORG_FILES)

BUILT_SOURCES = $(C_FILES) $(F_FILES) $(FH_FUNC_FILES) $(FH_TYPE_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES) $(qmckl_f) $(qmckl_h)
CLEANFILES += $(BUILT_SOURCES) $(C_TEST_FILES) $(F_TEST_FILES) $(TANGLED_FILES) $(C_TEST_FILES) $(F_TEST_FILES) $(qmckl_f) $(qmckl_h)


## Silent
## ======

tangle_verbose   = $(tangle_verbose_@AM_V@)
tangle_verbose_  = $(tangle_verbose_@AM_DEFAULT_V@)
tangle_verbose_0 = @echo "  TANGLE   $<";

cat_h_verbose   = $(cat_h_verbose_@AM_V@)
cat_h_verbose_  = $(cat_h_verbose_@AM_DEFAULT_V@)
cat_h_verbose_0 = @echo "  HEADER   $@";


## Rules
## =====

SUFFIXES = .f90 .h .org .c _f.f90 _func.h _type.h _private_func.h _private_type.h

$(qmckl_h): $(H_FUNC_FILES) $(H_TYPE_FILES)
	$(cat_h_verbose)$(srcdir)/tools/build_qmckl_h.sh

$(qmckl_f): $(FH_FUNC_FILES) $(FH_TYPE_FILES)
	$(cat_h_verbose)$(srcdir)/tools/build_qmckl_f.sh

# vim: syntax=automake noet

