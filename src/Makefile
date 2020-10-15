CC=gcc
CFLAGS=

FC=gfortran
FFLAGS=

ORG_SOURCE_FILES=qmckl_context.org
OBJECT_FILES=$(patsubst %.org,%.o,$(ORG_SOURCE_FILES))

.PHONY: clean

all: $(OBJECT_FILES)

%.c %.h: %.org
	emacs --quick --no-init-file --batch --eval "(require 'org)" --eval '(org-babel-tangle-file "$^")'

%.c %.h %_f.f90: %.org
	emacs --quick --no-init-file --batch --eval "(require 'org)" --eval '(org-babel-tangle-file "$^")'

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $*.o

clean:
	rm -f qmckl_*.f90 qmckl_*.c qmckl_*.o qmckl_*.h


