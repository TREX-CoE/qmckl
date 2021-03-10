#!/bin/bash 

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


# Create documentation
cd ${QMCKL_ROOT}/src

function extract_doc()
{
    HTML=${1%.org}.html 
    if [[ -f ${QMCKL_ROOT}/docs/$HTML && $1 -ot ${QMCKL_ROOT}/docs/$HTML ]]
    then return
    fi
    emacs --batch \
          --load ${QMCKL_ROOT}/docs/htmlize.el \
          --load ${QMCKL_ROOT}/tools/config_doc.el  \
          $i \
          --load ${QMCKL_ROOT}/tools/config_tangle.el  \
          -f org-html-export-to-html || break
    mv $HTML ${QMCKL_ROOT}/docs 
}

for i in *.org
do
echo
echo "=======  $i ======="
   extract_doc $i
done

if [[ $? -eq 0 ]]
then
  cd ${QMCKL_ROOT}/docs
  rm -f index.html
  ln README.html index.html
  exit 0
else
  exit 2
fi


