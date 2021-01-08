#!/bin/bash
INPUT=$1

if [[ -f ../docs/htmlize.el ]]
then
    emacs --batch --load ../docs/htmlize.el --load ../docs/config.el $INPUT -f org-html-export-to-html
else
    emacs --batch --load ../docs/config.el $INPUT -f org-html-export-to-html
fi

mv index.html ../docs


