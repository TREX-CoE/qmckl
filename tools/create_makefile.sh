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



# Then we create the list of ~*.o~ files to be created, for library
# functions:


OBJECTS=""
for i in $(ls qmckl_*.c) ; do
    FILE=${i%.c}
    OBJECTS="${OBJECTS} ${FILE}.o"
done >> $OUTPUT

for i in $(ls qmckl_*_f.f90) ; do
    FILE=${i%.f90}
    OBJECTS="${OBJECTS} ${FILE}.o"
done >> $OUTPUT



# for tests in C:


TESTS=""
for i in $(ls test_qmckl_*.c) ; do
    FILE=${i%.c}.o
    TESTS="${TESTS} ${FILE}"
done >> $OUTPUT



# and for tests in Fortran:


TESTS_F=""
for i in $(ls test_qmckl_*_f.f90) ; do
    FILE=${i%.f90}.o
    TESTS_F="${TESTS_F} ${FILE}"
done >> $OUTPUT



# Finally, we append the rules to the Makefile


cat << EOF > ${OUTPUT}
CC=$CC
CFLAGS=$CFLAGS -I../munit/ 

FC=$FC
FFLAGS=$FFLAGS
OBJECT_FILES=$OBJECTS
TESTS=$TESTS
TESTS_F=$TESTS_F

LIBS=$LIBS

libqmckl.so: \$(OBJECT_FILES)
	\$(CC) -shared \$(OBJECT_FILES) -o libqmckl.so
           
%.o: %.c 
	\$(CC) \$(CFLAGS) -c \$*.c -o \$*.o


qmckl_f.o: ../include/qmckl_f.f90
	\$(FC) \$(FFLAGS) -c ../include/qmckl_f.f90 -o qmckl_f.o

%.o: %.f90 qmckl_f.o
	\$(FC) \$(FFLAGS) -c \$*.f90 -o \$*.o


test_qmckl: test_qmckl.c libqmckl.so \$(TESTS) \$(TESTS_F)
	\$(CC) \$(CFLAGS) -Wl,-rpath,$PWD -L. \
	../munit/munit.c \$(TESTS) \$(TESTS_F) -lqmckl \$(LIBS) test_qmckl.c -o test_qmckl

test: test_qmckl
	./test_qmckl

.PHONY: test
EOF
