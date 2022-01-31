#!/bin/bash
#
# Installs the htmlize Emacs plugin

${abs_srcdir}/tools/missing git clone "https://github.com/TREX-CoE/emacs-htmlize"
mv emacs-htmlize/htmlize.el $1
rm -rf emacs-htmlize


