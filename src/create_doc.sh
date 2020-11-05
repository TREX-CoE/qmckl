#!/bin/bash
INPUT=$1
#emacs merged_qmckl.org --batch --eval "(require 'htmlize)" -f org-html-export-to-html --kill
emacs \
     $INPUT \
     --batch \
     --eval "(package-initialize)"  \
     -f org-html-export-to-html \
     --kill

mv *.html ../docs


