#!/bin/bash
#
# Installs the htmlize Emacs plugin

./missing git clone "https://github.com/hniksic/emacs-htmlize"
mv emacs-htmlize/htmlize.el $1
rm -rf emacs-htmlize


