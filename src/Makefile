CC=gcc
CFLAGS=-fPIC -fexceptions -Wall -Werror -Wpedantic -Wextra -g

FC=gfortran
FFLAGS=-fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wuninitialized  -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan


export CC CFLAGS FC FFLAGS

ORG_SOURCE_FILES=$(wildcard qmckl*.org) test_qmckl.org
OBJECT_FILES=$(filter-out $(EXCLUDED_OBJECTS), $(patsubst %.org,%.o,$(ORG_SOURCE_FILES)))

.PHONY: clean
.SECONDARY: # Needed to keep the produced C and Fortran files

libqmckl.so: Makefile.generated
	$(MAKE) -f Makefile.generated

test: Makefile.generated
	$(MAKE) -f Makefile.generated test

clean:
	rm -f qmckl.h test_qmckl_* qmckl_*.f90 qmckl_*.c qmckl_*.o qmckl_*.h Makefile.generated libqmckl.so

Makefile.generated: $(ORG_SOURCE_FILES) Makefile create_makefile.sh
	./create_makefile.sh $(ORG_SOURCE_FILES)

