LC_MESSAGES=C
export LC_MESSAGES

SUBDIRS =

libqmckl_la_version = $(VERSION_MAJOR).$(VERSION_MINOR)

#lib_LTLIBRARIES = libqmckl-$(libqmckl_so_version).la
lib_LTLIBRARIES = libqmckl.la

include generated.mk
qmckl_h = $(top_srcdir)/include/qmckl.h
qmckl_f = $(top_srcdir)/share/qmckl/fortran/qmckl_f.f90

BUILT_SOURCES= $(C_FILES) $(F_FILES) $(FH_FUNC_FILES) $(FH_TYPE_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES)

libqmckl_la_CPPFLAGS = -I$(srcdir) -I$(top_srcdir)/include
libqmckl_la_SOURCES = $(C_FILES) $(F_FILES) $(qmckl_h) $(qmckl_f) $(H_PRIVATE_FUNC_FILES) $(H_PRIVATE_TYPE_FILES)

tangle_verbose   = $(tangle_verbose_@AM_V@)
tangle_verbose_  = $(tangle_verbose_@AM_DEFAULT_V@)
tangle_verbose_0 = @echo "  TANGLE   $@";

cat_h_verbose   = $(cat_h_verbose_@AM_V@)
cat_h_verbose_  = $(cat_h_verbose_@AM_DEFAULT_V@)
cat_h_verbose_0 = @echo "  HEADER   $@";



export top_srcdir

Makefile.in: generated.mk

generated.mk: $(ORG_FILES)
	$(top_srcdir)/tools/create_makefile.sh

qmckl.h qmckl_f.f90: $(FH_FUNC_FILES) $(H_FUNC_FILES) $(H_TYPE_FILES)
	$(cat_h_verbose)$(top_srcdir)/tools/build_qmckl_h.sh

SUFFIXES = .f90 .h .org .c

.f90: $(ORG_FILES)
	$(tangle_verbose)$(top_srcdir)/tools/tangle.sh $(ORG_FILES)

.h: $(ORG_FILES)
	$(tangle_verbose)$(top_srcdir)/tools/tangle.sh $(ORG_FILES)

.org.c:
	$(tangle_verbose)$(top_srcdir)/tools/tangle.sh $<

$(qmckl_f): qmckl_f.f90
	mv qmckl_f.f90 $(qmckl_f)

$(qmckl_h): qmckl.h
	mv qmckl.h $(qmckl_h)

