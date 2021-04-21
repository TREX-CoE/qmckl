# Use POSIX-compliant Makefiles
.POSIX:

# Clear suffix list
.SUFFIXES:

package = qmckl
version = 0.1-alpha
tarname = $(package)
distdir = $(tarname)-$(version)
prefix  = /usr/local

QMCKL_ROOT=$(CURDIR)
shared_lib=$(QMCKL_ROOT)/lib/libqmckl.so
static_lib=$(QMCKL_ROOT)/lib/libqmckl.a
qmckl_h=$(QMCKL_ROOT)/include/qmckl.h
qmckl_f=$(QMCKL_ROOT)/share/qmckl/fortran/qmckl_f.f90

export prefix shared_lib static_lib qmckl_h qmckl_f


all clean doc check install uninstall:
	$(MAKE) -C src $@

dist: $(distdir).tar.gz


$(distdir).tar.gz: $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $@
	rm -rf $(distdir)


$(distdir): $(qmckl_h) $(qmckl_f) $(static_lib) $(shared_lib) src/Makefile.generated doc FORCE
	mkdir -p $(distdir)
	mkdir -p $(distdir)/munit
	mkdir -p $(distdir)/src
	mkdir -p $(distdir)/include
	mkdir -p $(distdir)/share/qmckl/fortran
	mkdir -p $(distdir)/share/qmckl/doc/html/
	mkdir -p $(distdir)/share/qmckl/doc/text/
	mkdir -p $(distdir)/man
	cp munit/munit.h munit/munit.c $(distdir)/munit/
	cp src/*.c src/*.h src/*.f90 $(distdir)/src/
	cp src/Makefile.generated $(distdir)/src/Makefile
	cp include/* $(distdir)/include
	cp Makefile $(distdir)/
	cp docs/*.html $(distdir)/share/qmckl/doc/html/
	cp docs/*.css  $(distdir)/share/qmckl/doc/html/
	cp docs/*.txt  $(distdir)/share/qmckl/doc/text/
	cp share/qmckl/fortran/*  $(distdir)/share/qmckl/fortran
	mkdir -p $(distdir)/lib


FORCE:
	rm -f -- $(distdir).tar.gz
	rm -rf -- $(distdir)


distcheck: $(distdir).tar.gz
	gzip -cd $(distdir).tar.gz | tar xvf -
	cd $(distdir) && $(MAKE) all check
	rm $(distdir)/lib/libqmckl.so $(distdir)/include/qmckl.h \
		$(distdir)/include/qmckl_f.f90
	cd $(distdir) && $(MAKE) clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz is ready for distribution."


$(qmckl_h) $(qmckl_f) $(static_lib) $(shared_lib) src/Makefile.generated:
	$(MAKE) -C src $@ 



.PHONY: all clean dist doc install uninstall FORCE
