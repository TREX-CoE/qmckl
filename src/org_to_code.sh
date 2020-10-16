#!/usr/bin/env bash

if [[ -z $1 ]] ; then
    echo "Usage: $0 <FILE.org>"
    exit 1;
fi

if [[ -z $6 ]] ; then
# Few file to tangle

    for INPUT in $@ ; do
	emacs \
	    --quick \
	    --no-init-file \
	    --batch \
	    --eval "(require 'org)" \
     	    --eval "(org-babel-tangle-file \"$INPUT\")" 
    done
	
else
# Multiple files to tangle, so we use the emacs server to speed up thing

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
fi
