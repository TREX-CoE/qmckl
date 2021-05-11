#!/bin/bash
# Extracts documentation from an org-mode file.

set -e

srcdir=$PWD

readonly DOCS=${srcdir}/share/doc/qmckl/
readonly ORG=${srcdir}/org/
readonly HTMLIZE=${DOCS}/html/htmlize.el
readonly CONFIG_DOC=${srcdir}/tools/config_doc.el
readonly CONFIG_TANGLE=${srcdir}/tools/config_tangle.el

function extract_doc()
{
  org=$1
  local_html=${org%.org}.html
  local_text=${org%.org}.txt
  html=${DOCS}/html/$(basename ${org%.org}.html)
  text=${DOCS}/text/$(basename ${org%.org}.txt)

  ./missing emacs --batch  \
        --load ${HTMLIZE}          \
        --load ${CONFIG_DOC}       \
        ${org}                     \
        --load ${CONFIG_TANGLE}    \
        -f org-html-export-to-html \
        -f org-ascii-export-to-ascii 

  mv ${local_html} ${html}
  mv ${local_text} ${text}
}

for i in $@
do
    exported=${i%.org}.exported
    NOW=$(date +"%m%d%H%M.%S")
    extract_doc ${i} &> $exported

    # Make log file older than the exported files
    touch -t ${NOW} $exported
done


