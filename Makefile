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

docdir=$(prefix)/share/doc/qmckl
libdir=$(prefix)/lib
includedir=$(prefix)/include
fortrandir=$(prefix)/share/qmckl/fortran

export prefix shared_lib static_lib qmckl_h qmckl_f


all clean doc install uninstall check:
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
	mkdir -p $(distdir)/share/doc/qmckl/html/
	mkdir -p $(distdir)/share/doc/qmckl/text/
	mkdir -p $(distdir)/man
	cp $(QMCKL_ROOT)/munit/munit.h munit/munit.c $(distdir)/munit/
	cp $(QMCKL_ROOT)/src/*.c src/*.h src/*.f90 $(distdir)/src/
	cp $(QMCKL_ROOT)/src/Makefile.generated $(distdir)/src/Makefile
	cp $(qmckl_h) $(distdir)/include
	cp $(QMCKL_ROOT)/Makefile $(distdir)/
	cp $(QMCKL_ROOT)/share/doc/qmckl/html/*.html $(distdir)/share/doc/qmckl/html/
	cp $(QMCKL_ROOT)/share/doc/qmckl/html/*.css  $(distdir)/share/doc/qmckl/html/
	cp $(QMCKL_ROOT)/share/doc/qmckl/text/*.txt  $(distdir)/share/doc/qmckl/text/
	cp $(qmckl_f) $(distdir)/share/qmckl/fortran/
	mkdir -p $(distdir)/lib


FORCE:
	rm -f -- $(distdir).tar.gz
	rm -rf -- $(distdir)


distcheck: $(distdir).tar.gz
	gzip -cd $(distdir).tar.gz | tar xvf -
	cd $(distdir) && $(MAKE) all 
	cd $(distdir) && $(MAKE) check
	cd $(distdir) && $(MAKE) prefix=$${PWD}/_inst install
	cd $(distdir) && $(MAKE) prefix=$${PWD}/_inst uninstall
	@remaining="`find $${PWD}/$(distdir)/_inst -type f | wc -l`" ;\
		if test "$${remaining}" -ne 0; then \
			echo "*** $${remaining} file(s) remaining in stage directory"; \
			exit 1; \
		fi
	cd $(distdir) && $(MAKE) clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz is ready for distribution."


$(qmckl_h) $(qmckl_f) $(static_lib) $(shared_lib):
	$(MAKE) -C src $@ 

src/Makefile.generated: 
	$(MAKE) -C src Makefile.generated

veryclean: FORCE clean


.PHONY: all veryclean clean dist doc install uninstall FORCE
