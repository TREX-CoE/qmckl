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



if VFC_CI
AM_LDFLAGS=-lvfc_probes -lvfc_probes_f
endif


ACLOCAL_AMFLAGS = -I m4

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

SUBDIRS =
CLEANFILES = qmckl.mod
EXTRA_DIST =

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pkgconfig/qmckl.pc

qmckl_h = $(srcdir)/include/qmckl.h
include_HEADERS = $(qmckl_h)

qmckl_f       = share/qmckl/fortran/qmckl_f.f90
test_qmckl_f  = tests/qmckl_f.f90
test_qmckl_fo = tests/qmckl_f.o
src_qmckl_f   = src/qmckl_f.f90
src_qmckl_fo  = src/qmckl_f.o
header_tests  = tests/chbrclf.h

fortrandir = $(datadir)/$(PACKAGE_NAME)/fortran/
dist_fortran_DATA = $(qmckl_f)

AM_CPPFLAGS = -I$(srcdir)/src -I$(srcdir)/include

lib_LTLIBRARIES = src/libqmckl.la
src_libqmckl_la_SOURCES = $(qmckl_h) $(src_qmckl_f) $(C_FILES) $(F_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES) $(header_tests)

export qmckl_f qmckl_h srcdir

CLEANFILES+=$(test_qmckl_f) $(src_qmckl_f) $(test_qmckl_o) $(src_qmckl_o)

htmlize_el=share/doc/qmckl/html/htmlize.el

include generated.mk


ln_s_verbose   = $(ln_s_verbose_@AM_V@)
ln_s_verbose_  = $(ln_s_verbose_@AM_DEFAULT_V@)
ln_s_verbose_0 = @echo "  LN_S     $<";

htmldir = $(docdir)/html
dist_html_DATA = $(HTML_FILES) \
        share/doc/qmckl/html/index.html \
        $(srcdir)/share/doc/qmckl/html/qmckl.css

textdir = $(docdir)/text
dist_text_DATA = $(TEXT_FILES)

$(test_qmckl_fo): $(test_qmckl_f)
$(test_qmckl_f): $(srcdir)/$(qmckl_f)
	cp $(srcdir)/$(qmckl_f) $(test_qmckl_f)

$(src_qmckl_fo): $(src_qmckl_f)
$(src_qmckl_f): $(srcdir)/$(qmckl_f)
	cp $(srcdir)/$(qmckl_f) $(src_qmckl_f)

share/doc/qmckl/html/index.html: share/doc/qmckl/html/README.html
	$(ln_s_verbose)cd share/doc/qmckl/html/ && \
        rm -rf index.html && \
	$(LN_S) README.html index.html


html: $(dist_html_DATA)
text: $(dist_text_DATA)
doc: html text



if QMCKL_DEVEL

dist_src_DATA = $(ORG_FILES) $(TANGLED_FILES) $(EXPORTED_FILES)

BUILT_SOURCES = $(C_FILES) $(F_FILES) $(FH_FUNC_FILES) $(FH_TYPE_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES) $(qmckl_f) $(qmckl_h) $(header_tests)

CLEANFILES += $(BUILT_SOURCES) $(C_TEST_FILES) $(F_TEST_FILES) $(TANGLED_FILES) $(C_TEST_FILES) $(F_TEST_FILES) $(qmckl_f) $(qmckl_h) $(HTML_FILES) $(TEXT_FILES) share/doc/qmckl/html/index.html $(EXPORTED_FILES) $(header_tests)

EXTRA_DIST +=                 \
    tools/build_doc.sh        \
    tools/build_makefile.py   \
    tools/build_qmckl_f.sh    \
    tools/build_qmckl_h.sh    \
    tools/config_doc.el       \
    tools/config_tangle.el    \
    tools/init.el             \
    tools/install_htmlize.sh  \
    tools/lib.org             \
    tools/tangle.sh           \
    tools/theme.setup


## Silent
## ======

export_verbose   = $(export_verbose_@AM_V@)
export_verbose_  = $(export_verbose_@AM_DEFAULT_V@)
export_verbose_0 = @echo "  DOC      $@";

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

$(htmlize_el):
	$(srcdir)/tools/install_htmlize.sh $(htmlize_el)

tests/chbrclf.h: $(qmckl_h)


generated.mk: $(ORG_FILES)
	$(PYTHON) $(srcdir)/tools/build_makefile.py

cppcheck: cppcheck.out

cppcheck.out: $(qmckl_h)
	cd src/ && \
	cppcheck --addon=cert -q --error-exitcode=0  \
          --enable=all --suppress="unusedStructMember"\
	  --suppress="unusedFunction" \
	  --suppress="missingIncludeSystem" \
          --language=c --std=c99 -rp --platform=unix64  \
          -I../include *.c *.h 2>../$@


.PHONY: cppcheck

endif
