CC=gcc
CFLAGS=

FC=gfortran
FFLAGS=

SOURCE_FILES=$(wildcard *.org)

%.c %_f.f90 %.h: %.c.org
	emacs --quick --no-init-file --batch --eval "(require 'org)" --eval '(org-babel-tangle-file "$^")'

%.o: %.c 
	$(CC) $(CFLAGS) -c $*.c -o $*.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $*.o


