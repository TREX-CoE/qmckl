#!/bin/bash 

INPUT=$1
SRC=$PWD


# Install htmlize if needed
[[ -f ../docs/htmlize.el ]] || (
    cd ../docs/
    git clone https://github.com/hniksic/emacs-htmlize 
    cp emacs-htmlize/htmlize.el .
    rm -rf  emacs-htmlize
    cd -
)

[[ -f ../docs/htmlize.el ]] || exit 1


# Switch to TMPDIR for easy cleanup
TMPDIR=$(mktemp -d)
./merge_org.sh $TMPDIR/$INPUT
cd $TMPDIR


# Create documentation
emacs --batch \
      --load ${SRC}/../docs/htmlize.el \
      --load ${SRC}/../toold/init.el  \
      $INPUT -f org-html-export-to-html

if [[ $? -eq 0 ]]
then
  rm -rf $TMPDIR
  exit 0
else
  mv index.html ${SRC}/../docs/
  rm -rf $TMPDIR
  exit 2
fi


