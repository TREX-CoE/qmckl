#!/bin/bash

OUTPUT=Makefile.generated

# Tangle org files

emacsclient -a "" \
            --socket-name=org_to_code \
            --eval "(require 'org)"

for INPUT in $@ ; do
    emacsclient \
	--no-wait \
	--socket-name=org_to_code \
     	--eval "(org-babel-tangle-file \"$INPUT\")"
done

emacsclient \
    --no-wait \
    --socket-name=org_to_code \
    --eval '(kill-emacs)'



# Create the list of *.o files to be created

OBJECTS=""
for i in $(ls qmckl_*.c) ; do
    FILE=${i%.c}
    OBJECTS="${OBJECTS} ${FILE}.o"
done >> $OUTPUT

for i in $(ls qmckl_*.f90) ; do
    FILE=${i%.f90}
    OBJECTS="${OBJECTS} ${FILE}.o"
done >> $OUTPUT

TESTS=""
for i in $(ls test_qmckl_*.c) ; do
    FILE=${i%.c}.o
    TESTS="${TESTS} ${FILE}"
done >> $OUTPUT

TESTS_F=""
for i in $(ls test_qmckl_*.f90) ; do
    FILE=${i%.f90}.o
    TESTS_F="${TESTS_F} ${FILE}"
done >> $OUTPUT


# Write the Makefile

cat << EOF > $OUTPUT
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

%.o: %.f90 
	\$(FC) \$(FFLAGS) -c \$*.f90 -o \$*.o

test_qmckl: test_qmckl.c libqmckl.so \$(TESTS) \$(TESTS_F)
	\$(CC) \$(CFLAGS) -Wl,-rpath,$PWD -L. \
	../munit/munit.c \$(TESTS) \$(TESTS_F) -lqmckl \$(LIBS) test_qmckl.c -o test_qmckl

test: test_qmckl
	./test_qmckl

.PHONY: test
EOF

for i in $(ls qmckl_*.c) ; do
    FILE=${i%.c}
    echo "${FILE}.o: ${FILE}.c " *.h
done >> $OUTPUT

for i in $(ls qmckl_*.f90) ; do
    FILE=${i%.f90}
    echo "${FILE}.o: ${FILE}.f90"
done >> $OUTPUT

for i in $(ls test_qmckl_*.c) ; do
    FILE=${i%.c}
    echo "${FILE}.o: ${FILE}.c qmckl.h" 
done >> $OUTPUT


for i in $(ls test_qmckl*.f90) ; do
    FILE=${i%.f90}
    echo "${FILE}.o: ${FILE}.f90"
done >> $OUTPUT


