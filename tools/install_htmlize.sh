#!/bin/bash
#
# Installs the htmlize Emacs plugin

readonly HTMLIZE_EL=$1

# If the file already present - exit the script
if [ -f ${HTMLIZE_EL} ]
then
	exit 0	
fi

DOC_ROOT=${HTMLIZE_EL%/*}
readonly EMACS_HTMLIZE=${DOC_ROOT}/emacs-htmlize/htmlize.el

# Case 1: submodule cloned but htmlize.el is absent in DOC_ROOT
if [ -f ${EMACS_HTMLIZE} ]
then
	cp ${EMACS_HTMLIZE} ${HTMLIZE_EL}
else
# Case 2: submodule has not been cloned
	${abs_srcdir}/tools/missing git clone "https://github.com/TREX-CoE/emacs-htmlize"
	mv emacs-htmlize/htmlize.el ${HTMLIZE_EL}
	rm -rf emacs-htmlize
fi

