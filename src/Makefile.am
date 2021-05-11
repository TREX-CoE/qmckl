LC_MESSAGES=C
export LC_MESSAGES

SUBDIRS =

libqmckl_la_version = $(VERSION_MAJOR).$(VERSION_MINOR)

#lib_LTLIBRARIES = libqmckl-$(libqmckl_so_version).la
lib_LTLIBRARIES = libqmckl.la


include generated.mk

BUILT_SOURCES= $(C_FILES) $(F_FILES) $(FH_FUNC_FILES) $(FH_TYPE_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES)

libqmckl_la_CPPFLAGS = -I$(srcdir) -I$(top_srcdir)/include
libqmckl_la_SOURCES = qmckl.h qmckl_f.f90 $(C_FILES) $(F_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES)

dist_src_DATA = $(ORG_FILES)

tangle_verbose   = $(tangle_verbose_@AM_V@)
tangle_verbose_  = $(tangle_verbose_@AM_DEFAULT_V@)
tangle_verbose_0 = @echo "  TANGLE   $@";

cat_h_verbose   = $(cat_h_verbose_@AM_V@)
cat_h_verbose_  = $(cat_h_verbose_@AM_DEFAULT_V@)
cat_h_verbose_0 = @echo "  HEADER   $@";

SUFFIXES = .f90 .h .org .c _f.f90 _func.h _type.h _private_func.h _private_type.h

export top_srcdir

if QMCKL_DEVEL

Makefile.in: generated.mk

generated.mk: $(ORG_FILES)
	cd $(top_srcdir) && ./tools/create_makefile.sh

qmckl.h qmckl_f.f90: $(FH_FUNC_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES)
	$(cat_h_verbose)$(top_srcdir)/tools/build_qmckl_h.sh

.org_private_func.h .org_private_type.h .org_type.h .org_func.h .org.h .org_f.f90 .org.c:
	$(tangle_verbose)$(top_srcdir)/tools/tangle.sh $<

.h .f90 _f.f90 _func.h _type.h _private_func.h _private_type.h .c: $(ORG_FILES)

CLEANFILES = $(BUILT_SOURCES) $(C_TEST_FILES) $(F_TEST_FILES) qmckl_f.f90 qmckl.h
endif
