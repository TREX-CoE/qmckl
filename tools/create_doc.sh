#!/bin/bash 

INPUT=merged.org
if [[ -z $QMCKL_ROOT ]]
then
  print "QMCKL_ROOT is not defined"
  exit 1
fi


# Install htmlize if needed
[[ -f ${QMCKL_ROOT}/docs/htmlize.el ]] || (
    cd ${QMCKL_ROOT}/docs/
    git clone https://github.com/hniksic/emacs-htmlize 
    cp emacs-htmlize/htmlize.el .
    rm -rf emacs-htmlize
    cd -
)

[[ -f ${QMCKL_ROOT}/docs/htmlize.el ]] || exit 1


# Switch to TMPDIR for easy cleanup
TMPDIR=$(mktemp -d)
${QMCKL_ROOT}/tools/merge_org.sh $TMPDIR/$INPUT
cd $TMPDIR


# Create documentation
emacs --batch \
      --load ${QMCKL_ROOT}/docs/htmlize.el \
      --load ${QMCKL_ROOT}/docs/config.el  \
      $INPUT -f org-html-export-to-html

if [[ $? -eq 0 ]]
then
  mv index.html ${QMCKL_ROOT}/docs/
  rm -rf $TMPDIR
  exit 0
else
  rm -rf $TMPDIR
  exit 2
fi


