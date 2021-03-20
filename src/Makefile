# Header                                                          :noexport:


# This file was created by tools/Building.org

# Dependencies


LIBS=-lpthread

# Variables


QMCKL_ROOT=$(shell dirname $(CURDIR))

export CC CFLAGS FC FFLAGS LIBS QMCKL_ROOT

ORG_SOURCE_FILES=$(wildcard *.org)
OBJECT_FILES=$(filter-out $(EXCLUDED_OBJECTS), $(patsubst %.org,%.o,$(ORG_SOURCE_FILES)))
INCLUDE=-I$(QMCKL_ROOT)/include/

# Compiler options

#  GNU, Intel and LLVM compilers are supported. Choose here:


COMPILER=GNU
#COMPILER=INTEL
#COMPILER=LLVM

# GNU


ifeq ($(COMPILER),GNU)
#----------------------------------------------------------
CC=gcc -g
CFLAGS=-fPIC  $(INCLUDE) \
	-fexceptions -Wall -Werror -Wpedantic -Wextra -fmax-errors=3

FC=gfortran -g
FFLAGS=-fPIC  $(INCLUDE) \
	-fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising \
	-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation \
	-Wreal-q-constant -Wuninitialized  -fbacktrace -finit-real=nan \
	-ffpe-trap=zero,overflow,underflow

LIBS+=-lgfortran -lm
#----------------------------------------------------------
endif

# Intel


ifeq ($(COMPILER),INTEL)
#----------------------------------------------------------
CC=icc -xHost
CFLAGS=-fPIC -g -O2 $(INCLUDE)

FC=ifort -xHost
FFLAGS=-fPIC -g -O2 $(INCLUDE)

LIBS+=-lm -lifcore -lirc
#----------------------------------------------------------
CC=icc -xHost
endif

# LLVM


ifeq ($(COMPILER),LLVM)
#----------------------------------------------------------
CC=clang
CFLAGS=-fPIC -g -O2 $(INCLUDE)

FC=flang
FFLAGS=fPIC -g -O2 $(INCLUDE)

LIBS+=-lm
#----------------------------------------------------------
endif

# Rules

#    The source files are created during the generation of the file ~Makefile.generated~.


.PHONY: clean
.SECONDARY: # Needed to keep the produced C and Fortran files

libqmckl.so: ../include/qmckl.h
	$(MAKE) -f Makefile.generated

../include/qmckl.h: Makefile.generated
	../tools/build_qmckl_h.sh

test: libqmckl.so
	$(MAKE) -f Makefile.generated test

doc: $(ORG_SOURCE_FILES)
	$(QMCKL_ROOT)/tools/build_doc.sh

clean:
	$(RM)	test_qmckl_* test_qmckl.c test_qmckl \
		qmckl_*.f90 qmckl_*.c qmckl_*.o qmckl_*.h \
		Makefile.generated libqmckl.so *.html *.mod

Makefile.generated: Makefile $(QMCKL_ROOT)/tools/create_makefile.sh  $(ORG_SOURCE_FILES)
	$(QMCKL_ROOT)/tools/create_makefile.sh
