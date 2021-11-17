#!/bin/bash
# Extracts documentation from an org-mode file.

# Check required environment variables
# ------------------------------------

if [[ -z ${srcdir} ]] ; then
   echo "Error: srcdir environment variable is not defined"
   exit 1
fi

if [[ -z ${top_builddir} ]] ; then
   echo "Error: srcdir environment variable is not defined"
   exit 1
fi

readonly DOCS=${top_builddir}/share/doc/qmckl/
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

  ${srcdir}/tools/missing emacs --batch  \
        --load="${HTMLIZE}"           \
        --load="${CONFIG_DOC}"        \
        ${org}                           \
        --load="${CONFIG_TANGLE}"     \
        -f org-html-export-to-html       \
        -f org-ascii-export-to-ascii 

  rc=$?
  if [ $rc -ne 0 ] ; then
    return $rc
  fi
  mv ${local_html} ${html}
  mv ${local_text} ${text}
}

for i in $@
do
    exported=${i%.org}.exported
    exported=${top_builddir}/src/$(basename $exported)
    NOW=$(date +"%m%d%H%M.%S")
    extract_doc ${i} &> $exported 
    rc=$?

    # Make log file older than the exported files
    touch -t ${NOW} $exported

    # Fail if tangling failed
    if [ $rc -ne 0 ] ; then
       cat $exported
       exit $rc
    fi
done


