CC=gcc -g
CFLAGS=-fPIC -fexceptions -Wall -Werror -Wpedantic -Wextra

FC=gfortran -g
FFLAGS=-fPIC  -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wuninitialized  -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan

LIBS=-lgfortran -lm
#
#CC=icc -xHost
#CFLAGS=-fPIC -g -O2
#
#FC=ifort -xHost
#FFLAGS=-fPIC -g -O2
#
#LIBS=-lm -lifcore -lirc


export CC CFLAGS FC FFLAGS LIBS

MERGED_ORG=merged_qmckl.org
ORG_SOURCE_FILES=$(wildcard *.org)
OBJECT_FILES=$(filter-out $(EXCLUDED_OBJECTS), $(patsubst %.org,%.o,$(ORG_SOURCE_FILES)))

.PHONY: clean
.SECONDARY: # Needed to keep the produced C and Fortran files

libqmckl.so: Makefile.generated
	$(MAKE) -f Makefile.generated

test: Makefile.generated
	$(MAKE) -f Makefile.generated test

$(MERGED_ORG): $(ORG_SOURCE_FILES)
	./merge_org.sh

doc:$(MERGED_ORG)
	./create_doc.sh $(MERGED_ORG)

clean:
	rm -f qmckl.h test_qmckl_* test_qmckl.c test_qmckl qmckl_*.f90 qmckl_*.c qmckl_*.o qmckl_*.h Makefile.generated libqmckl.so *.html *.fh *.mod

Makefile.generated: $(MERGED_ORG) Makefile create_makefile.sh
	./create_makefile.sh $(MERGED_ORG) 

