CC=gcc
CFLAGS=

SOURCE_FILES=$(wildcard *.org)

%.o: %.c.org
	- emacs --quick --no-init-file --batch --eval "(require 'org)" --eval '(org-babel-tangle-file "$^")'
	- $(CC) $(CFLAGS) -c $*.c -o $@
	- rm -f $*.c $*.h


