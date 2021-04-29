# Header                                                          :noexport:


# This file was created by tools/Building.org

.POSIX:

# Dependencies


LIBS=-lpthread

# Variables


QMCKL_ROOT=$(shell dirname $(CURDIR))

shared_lib=$(QMCKL_ROOT)/lib/libqmckl.so
static_lib=$(QMCKL_ROOT)/lib/libqmckl.a
qmckl_h=$(QMCKL_ROOT)/include/qmckl.h
qmckl_f=$(QMCKL_ROOT)/share/qmckl/fortran/qmckl_f.f90

export CC CFLAGS FC FFLAGS LIBS QMCKL_ROOT

ORG_SOURCE_FILES=$(wildcard *.org)
C_SOURCE_FILES=$(patsubst %.org,%.c,$(ORG_SOURCE_FILES))
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
#    The Makefile.generated is the one that will be distributed with the library.


.PHONY: clean shared static doc all check install uninstall
.SECONDARY: # Needed to keep the produced C and Fortran files

$(shared_lib) $(static_lib): $(qmckl_h) $(qmckl_f) Makefile.generated
	$(MAKE) -f Makefile.generated $@

install uninstall: Makefile.generated
	$(MAKE) -f Makefile.generated $@

$(qmckl_f) $(qmckl_h): Makefile.generated
	$(QMCKL_ROOT)/tools/build_qmckl_h.sh

shared: $(shared_lib)
static: $(static_lib)
all: shared static doc check

check: $(static_lib)
	$(MAKE) -f Makefile.generated check

cppcheck:
	cppcheck \
	  --addon=cert \
	  --enable=warning,style,performance,portability,information \
	  qmckl_*.c

doc: $(ORG_SOURCE_FILES)
	$(QMCKL_ROOT)/tools/build_doc.sh

clean:
	- $(MAKE) -f Makefile.generated clean
	- $(RM)	test_qmckl_* test_qmckl.c \
		$(qmckl_h) $(qmckl_f) \
		qmckl_*.f90 qmckl_*.c qmckl_*.h \
		Makefile.generated *.html *.txt

veryclean: clean FORCE
	- $(RM)	$(QMCKL_ROOT)/share/doc/qmckl/html/*.html \
		$(QMCKL_ROOT)/share/doc/qmckl/text/*.txt

Makefile.generated: Makefile $(QMCKL_ROOT)/tools/create_makefile.sh  $(ORG_SOURCE_FILES) ../tools/Building.org
	$(QMCKL_ROOT)/tools/create_makefile.sh


.SUFFIXES: .org .c

.org.c:
	$(QMCKL_ROOT)/tools/tangle.sh $<
