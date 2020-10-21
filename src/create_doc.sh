#!/bin/bash

# Tangle org files

emacsclient -a "" \
            --socket-name=org_to_code \
            --eval "(require 'org)"

for INPUT in $@ ; do
    echo $INPUT
    emacsclient \
	--no-wait \
	--socket-name=org_to_code \
     	--eval "(find-file \"$INPUT\")" \
     	--eval "(org-html-export-to-html)"
done

emacsclient \
    --no-wait \
    --socket-name=org_to_code \
    --eval '(kill-emacs)'

