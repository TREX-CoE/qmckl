#!/bin/bash
# Script to generate auto-generated Makefile
#   :PROPERTIES:
#   :header-args: :tangle create_makefile.sh :noweb  yes :shebang #!/bin/bash :comments org
#   :END:

#   This script generates the Makefile that compiles the library.
#   The ~OUTPUT~ variable contains the name of the generated Makefile,typically
#   =Makefile.generated=.


# This file was created by tools/Building.org



OUTPUT=Makefile.generated



# We start by tangling all the org-mode files.


${QMCKL_ROOT}/tools/tangle.sh *.org
../tools/build_qmckl_h.sh



# Then we create the list of ~*.o~ files to be created, for library
# functions:


OBJECTS="qmckl_f.o"
for i in $(ls qmckl_*.c qmckl_*f.f90) ; do
    FILE=${i%.*}
    OBJECTS+=" ${FILE}.o"
done >> $OUTPUT



# for tests in C:


TESTS=""
for i in $(ls test_qmckl_*.c) ; do
    FILE=${i%.c}
    TESTS+=" ${FILE}.o"
done >> $OUTPUT



# and for tests in Fortran:


TESTS_F=""
for i in $(ls test_qmckl_*_f.f90) ; do
    FILE=${i%.f90}
    TESTS_F+=" ${FILE}.o"
done >> $OUTPUT



# Finally, we append the rules to the Makefile


cat << EOF > ${OUTPUT}
.POSIX:
.SUFFIXES:

PREFIX=/usr/local

CC=$CC
CFLAGS=$CFLAGS -I../munit/

FC=$FC
FFLAGS=$FFLAGS

OBJECT_FILES=$OBJECTS
TESTS=$TESTS
TESTS_F=$TESTS_F

LIBS=$LIBS

QMCKL_ROOT=\$(shell dirname \$(CURDIR))
shared_lib=\$(QMCKL_ROOT)/lib/libqmckl.so
static_lib=\$(QMCKL_ROOT)/lib/libqmckl.a
qmckl_h=\$(QMCKL_ROOT)/include/qmckl.h
qmckl_f=\$(QMCKL_ROOT)/include/qmckl_f.f90
munit=\$(QMCKL_ROOT)/munit/munit.c 

shared: \$(shared_lib)
static: \$(static_lib)
all: shared static

\$(shared_lib): \$(OBJECT_FILES)
	\$(CC) -shared \$(OBJECT_FILES) -o \$(shared_lib)

\$(static_lib): \$(OBJECT_FILES)
	\$(AR) rcs \$(static_lib) \$(OBJECT_FILES)


# Test
  
qmckl_f.o: \$(qmckl_f)
	\$(FC) \$(FFLAGS) -c \$(qmckl_f) -o \$@

test_qmckl: test_qmckl.c \$(qmckl_h) \$(static_lib) \$(TESTS) \$(TESTS_F)
	\$(CC) \$(CFLAGS) \
	\$(munit) \$(TESTS) \$(TESTS_F) \$(static_lib) \$(LIBS) test_qmckl.c -o \$@

test_qmckl_shared: test_qmckl.c \$(qmckl_h) \$(shared_lib) \$(TESTS) \$(TESTS_F)
	\$(CC) \$(CFLAGS) -Wl,-rpath,$PWD/../lib -L../lib \
	\$(munit) \$(TESTS) \$(TESTS_F) -lqmckl \$(LIBS) test_qmckl.c -o \$@

check: test_qmckl test_qmckl_shared
	./test_qmckl

clean:
	\$(RM) -- *.o *.mod \$(shared_lib) \$(static_lib) test_qmckl

.SUFFIXES: .c .f90 .o

.c.o:
	\$(CC) \$(CFLAGS) -c \$*.c -o \$*.o

.f90.o: qmckl_f.o
	\$(FC) \$(FFLAGS) -c \$*.f90 -o \$*.o

.PHONY: check clean all
EOF
