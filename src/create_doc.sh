#!/bin/bash

# Tangle org files

emacsclient -a "" \
            --socket-name=org_to_code  \
            --eval "(load-file \"config.el\")"

for INPUT in $@ ; do
    echo $INPUT
    emacsclient \
        --no-wait \
        --socket-name=org_to_code \
        --eval "(find-file \"$INPUT\")" \
        --eval "(org-html-export-to-html)"
done
mv *.html ../docs

emacsclient \
    --no-wait \
    --socket-name=org_to_code \
    --eval '(kill-emacs)'

