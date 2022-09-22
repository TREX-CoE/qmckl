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
CLEANFILES = qmckl.mod qmckl_verificarlo_f.mod
EXTRA_DIST =

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pkgconfig/qmckl.pc

qmckl_h = include/qmckl.h
include_HEADERS = $(qmckl_h)

test_qmckl_f  = tests/qmckl_f.F90
test_qmckl_fo = tests/qmckl_f.o
src_qmckl_f   = src/qmckl_f.F90
src_qmckl_fo  = src/qmckl_f.o
header_tests  = tests/chbrclf.h tests/n2.h

fortrandir = $(datadir)/qmckl/fortran
fortran_DATA = $(src_qmckl_f)

QMCKL_TEST_DIR = $(abs_srcdir)/share/qmckl/test_data/

AM_CPPFLAGS  = -I$(top_builddir)/src -I$(top_builddir)/include
AM_CPPFLAGS += -I$(srcdir)/src -I$(srcdir)/include -I$(QMCKL_TEST_DIR)
AM_CPPFLAGS += -DQMCKL_TEST_DIR=\"$(QMCKL_TEST_DIR)\"


lib_LTLIBRARIES = src/libqmckl.la
src_libqmckl_la_SOURCES = $(qmckl_h) $(src_qmckl_f) $(C_FILES) $(F_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES)
src_libqmckl_la_LDFLAGS = $(LDFLAGS)

CLEANFILES+=$(test_qmckl_fo) $(src_qmckl_fo) $(test_qmckl_o) $(src_qmckl_o) $(FH_TYPE_FILES) $(FH_FUNC_FILES)


include generated.mk


ln_s_verbose   = $(ln_s_verbose_@AM_V@)
ln_s_verbose_  = $(ln_s_verbose_@AM_DEFAULT_V@)
ln_s_verbose_0 = @echo "  LN_S     $<";

htmldir = $(docdir)/html
textdir = $(docdir)/text
htmldir_loc = share/doc/qmckl/html
textdir_loc = share/doc/qmckl/text
htmlize_el = $(htmldir_loc)/htmlize.el
dist_html_DATA = $(HTML_FILES) \
        $(htmldir_loc)/qmckl.css \
        $(htmldir_loc)/index.html

dist_text_DATA = $(TEXT_FILES)

# Include Guix manifest in the source code distribution tarball
tools_qmckl_scm = tools/qmckl.scm
EXTRA_DIST += $(tools_qmckl_scm)

html-local: $(htmlize_el) $(dist_html_DATA)
text: $(htmlize_el) $(dist_text_DATA)
doc: html text

setup_py          =  $(srcdir)/python/setup.py
process_header_py =  $(srcdir)/python/src/process_header.py
test_py           =  $(srcdir)/python/test/test_api.py
qmckl_i           =  $(srcdir)/python/src/qmckl.i
numpy_i           =  $(srcdir)/python/src/numpy.i
qmckl_wrap_c    =  python/src/qmckl_wrap.c
qmckl_include_i =  python/src/qmckl_include.i
qmckl_py        =  python/src/qmckl.py

dist_python_DATA = $(setup_py) $(qmckl_py) $(qmckl_wrap_c) \
                   $(srcdir)/python/pyproject.toml \
                   $(srcdir)/python/requirements.txt \
                   $(srcdir)/python/README.md

python-install: $(qmckl_h) $(lib_LTLIBRARIES) $(dist_python_DATA)
	$(MKDIR_P) python/src
	cd python ; \
	[[ ! -f  pyproject.toml ]] && \
	    cp $(abs_srcdir)/python/{pyproject.toml,requirements.txt,README.md,setup.py} . ; \
	cp src/qmckl.py . ; \
	export QMCKL_INCLUDEDIR="$(prefix)/include" ; \
	export QMCKL_LIBDIR="$(prefix)/lib" ; \
	python3 -m pip install .

python-test: $(test_py)
	cd $(abs_srcdir)/python/test/ && \
		python test_api.py

.PHONY: python-test python-install cppcheck


if QMCKL_DEVEL

if VFC_CI
AM_LDFLAGS=-lvfc_probes -lvfc_probes_f
endif

EXTRA_DIST += $(ORG_FILES) $(TANGLED_FILES) $(EXPORTED_FILES) $(test_qmckl_f)

BUILT_SOURCES = $(C_FILES) $(F_FILES) $(FH_FUNC_FILES) $(FH_TYPE_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES) $(src_qmckl_f) $(test_qmckl_f) $(qmckl_h) $(header_tests) $(htmlize_el)

CLEANFILES += $(BUILT_SOURCES) $(C_TEST_FILES) $(F_TEST_FILES) $(TANGLED_FILES) $(C_TEST_FILES) $(F_TEST_FILES) $(src_qmckl_f) $(test_qmckl_f) $(qmckl_h) $(HTML_FILES) $(TEXT_FILES) $(EXPORTED_FILES) $(header_tests) $(htmlize_el)

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
tangle_verbose_0 = @echo "  ORG      $<";

cat_h_verbose   = $(cat_h_verbose_@AM_V@)
cat_h_verbose_  = $(cat_h_verbose_@AM_DEFAULT_V@)
cat_h_verbose_0 = @echo "  HEADER   $@";


## Rules
## =====

SUFFIXES = .F90 .h .org .c _f.F90 _func.h _type.h _private_func.h _private_type.h

$(test_qmckl_f): $(src_qmckl_f)
	cp $(src_qmckl_f) $(test_qmckl_f)

$(qmckl_h): $(H_FUNC_FILES) $(H_TYPE_FILES)
	$(cat_h_verbose)top_builddir=$(abs_top_builddir) srcdir=$(abs_srcdir) qmckl_h=$(qmckl_h) $(srcdir)/tools/build_qmckl_h.sh

$(src_qmckl_f): $(FH_FUNC_FILES) $(FH_TYPE_FILES)
	echo $(FH_TYPE_FILES) $(FH_FUNC_FILES)
	$(cat_h_verbose)top_builddir=$(abs_top_builddir) srcdir=$(abs_srcdir) src_qmckl_f=$(src_qmckl_f) $(srcdir)/tools/build_qmckl_f.sh

$(htmlize_el):
	$(MKDIR_P) $(htmldir_loc) $(textdir_loc)
	abs_srcdir=$(abs_srcdir) $(srcdir)/tools/install_htmlize.sh $(htmlize_el)

tests/chbrclf.h: $(qmckl_h)

tests/n2.h: $(qmckl_h)

cppcheck: cppcheck.out

cppcheck.out: $(qmckl_h)
	cd src/ && \
	cppcheck --addon=cert -q --error-exitcode=0  \
          --enable=all --suppress="unusedStructMember"\
	  --suppress="unusedFunction" \
	  --suppress="missingIncludeSystem" \
          --language=c --std=c99 -rp --platform=unix64  \
          -I$(srcdir)/include -I$(top_builddir)/include *.c *.h 2>../$@

$(qmckl_include_i):  $(qmckl_h) $(process_header_py)
	$(MKDIR_P) python/src
	python3 $(process_header_py) $(qmckl_h)
	mv qmckl_include.i $(qmckl_include_i)


$(qmckl_py): $(qmckl_i) $(qmckl_include_i)
	$(srcdir)/tools/missing swig -Iinclude -Ipython/src -python -py3 -builtin -o $(qmckl_wrap_c) $(qmckl_i) || :

$(qmckl_wrap_c): $(qmckl_py)


CLEANFILES += $(qmckl_wrap_c) \
	      $(qmckl_include_i) \
	      $(qmckl_py)

endif
