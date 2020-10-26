CC=gcc
CFLAGS=-fPIC -fexceptions -Wall -Werror -Wpedantic -Wextra -g

FC=gfortran
FFLAGS=-fPIC -g -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wuninitialized  -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan

LIBS=-lgfortran -lm

#CC=icc
#CFLAGS=-fPIC -g
#
#FC=ifort
#FFLAGS=-fPIC -g
#
#LIBS=-lm -lifcore -lirc 


export CC CFLAGS FC FFLAGS LIBS

ORG_SOURCE_FILES=$(wildcard qmckl*.org) test_qmckl.org
OBJECT_FILES=$(filter-out $(EXCLUDED_OBJECTS), $(patsubst %.org,%.o,$(ORG_SOURCE_FILES)))

.PHONY: clean
.SECONDARY: # Needed to keep the produced C and Fortran files

libqmckl.so: Makefile.generated
	$(MAKE) -f Makefile.generated

test: Makefile.generated
	$(MAKE) -f Makefile.generated test

doc:$(ORG_SOURCE_FILES)
	./create_doc.sh README.org $(ORG_SOURCE_FILES)

clean:
	rm -f qmckl.h test_qmckl_* test_qmckl.c qmckl_*.f90 qmckl_*.c qmckl_*.o qmckl_*.h Makefile.generated libqmckl.so *.html *.fh

Makefile.generated: $(ORG_SOURCE_FILES) Makefile create_makefile.sh
	./create_makefile.sh $(ORG_SOURCE_FILES)

