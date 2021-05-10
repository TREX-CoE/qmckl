#!/bin/bash
#
# This script tangles all the org-mode files in the src directory of QMCkl.
# It needs to be run from in the src directory.  It uses the config_tangle.el
# Emacs configuration file, which contains information required to compute the
# current file names using for example ~(eval c)~ to get the name of the
# produced C file. The org-mode file is not tangled if the last modification
# date of the org file is older than one of the tangled files.
# The =missing= script is used to check if emacs is present on the system.

if [[ $(basename $PWD) != "src" ]] ; then
        print "Error: $0 needs to be run from src directory"
        exit 1
fi


function tangle()
{
    local org_file=$1
    local c_file=${org_file%.org}.c
    local f_file=${org_file%.org}.f90

    if [[ ${org_file} -ot ${c_file} ]] ; then
        return
    elif [[ ${org_file} -ot ${f_file} ]] ; then
        return
    fi
    ../missing \
        emacs --batch ${org_file} \
         --load=${top_srcdir}/tools/config_tangle.el \
        -f org-babel-tangle
}

for i in $@
do
#    echo "--- ${i} ----"
    tangle ${i}
done
