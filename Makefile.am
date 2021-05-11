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


ACLOCAL_AMFLAGS = -I m4

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

SUBDIRS =
CLEANFILES =

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(srcdir)/qmckl.pc

qmckl_h = $(srcdir)/include/qmckl.h
include_HEADERS = $(qmckl_h)

qmckl_f  = $(srcdir)/share/qmckl/fortran/qmckl_f.f90
qmckl_fo = $(srcdir)/share/qmckl/fortran/qmckl_f.$(OBJEXT)
test_qmckl_f = $(srcdir)/tests/qmckl_f.f90
src_qmckl_f = $(srcdir)/src/qmckl_f.f90

fortrandir = $(datadir)/$(PACKAGE_NAME)/fortran/
fortran_DATA = $(qmckl_f)

AM_CPPFLAGS = -I$(srcdir)/src -I$(srcdir)/include
EXTRA_DIST = share/doc/qmckl/

lib_LTLIBRARIES = src/libqmckl.la
src_libqmckl_la_SOURCES = $(qmckl_h) $(src_qmckl_f) $(C_FILES) $(F_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES)

export srcdir qmckl_f qmckl_h

CLEANFILES+=$(test_qmckl_f)

include $(srcdir)/generated.mk

if QMCKL_DEVEL
include Makefile.dev.mk
else
include Makefile.dist.mk
endif

$(srcdir)/generated.mk: $(ORG_FILES)
	$(srcdir)/tools/build_makefile.sh

$(test_qmckl_f): $(qmckl_f)
	cp $(qmckl_f) $(test_qmckl_f)

$(src_qmckl_f): $(qmckl_f)
	cp $(qmckl_f) $(src_qmckl_f)

rpms srcrpm:
	for dir in $(SUBDIRS); do \
	  (cd $$dir && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1; \
	done

.PHONY: rpms srcrpm


