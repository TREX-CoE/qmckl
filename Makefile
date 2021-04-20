# Use POSIX-compliant Makefiles
.POSIX:

# Clear suffix list
.SUFFIXES:

package = qmckl
version = 0.1-alpha
tarname = $(package)
distdir = $(tarname)-$(version)


all clean check:
	$(MAKE) -C src $@

dist: $(distdir).tar.gz


$(distdir).tar.gz: $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $@
	rm -rf $(distdir)


$(distdir): include/qmckl.h include/qmckl_f.f90 src/Makefile.generated FORCE
	mkdir -p $(distdir)
	mkdir -p $(distdir)/munit
	mkdir -p $(distdir)/src
	mkdir -p $(distdir)/include
	cp munit/munit.h munit/munit.c $(distdir)/munit
	cp src/*.c src/*.h src/*.f90 $(distdir)/src
	cp src/Makefile.generated $(distdir)/src/Makefile
	cp include/* $(distdir)/include
	cp Makefile $(distdir)/
	mkdir -p $(distdir)/lib


FORCE:
	- rm -- $(distdir).tar.gz >/dev/null 2>&1
	- rm -rf -- $(distdir) >/dev/null 2>&1


distcheck: $(distdir).tar.gz
	gzip -cd $(distdir).tar.gz | tar xvf -
	cd $(distdir) && $(MAKE) all check
	rm $(distdir)/lib/libqmckl.so $(distdir)/include/qmckl.h \
		$(distdir)/include/qmckl_f.f90
	cd $(distdir) && $(MAKE) clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz is ready for distribution."


include/qmckl.h include/qmckl_f.f90 src/Makefile.generated:
	$(MAKE) -C src



.PHONY: all clean dist FORCE
