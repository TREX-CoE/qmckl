#!/bin/bash
INPUT=$1

emacs --batch --load ../docs/config.el $INPUT -f org-html-export-to-html

mv index.html ../docs


