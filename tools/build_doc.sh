#!/bin/bash
# Script to build the documentation.

set -e

srcdir=$PWD

readonly DOCS=${srcdir}/share/doc/qmckl/
readonly ORG=${srcdir}/org/
readonly HTMLIZE=${DOCS}/html/htmlize.el
readonly CONFIG_DOC=${srcdir}/tools/config_doc.el
readonly CONFIG_TANGLE=${srcdir}/tools/config_tangle.el

# Checks that all the defined global variables correspond to files.

for dir in ${DOCS}/html ${DOCS}/text ${ORG}
do
    if [[ ! -d ${dir} ]]
    then
        print "${dir} not found"
        exit 2
    fi
done

for file in ${CONFIG_DOC} ${CONFIG_TANGLE}
do
    if [[ ! -f ${file} ]]
    then
        print "${file} not found"
        exit 3
    fi
done




function install_htmlize()
{
#   Installs the htmlize Emacs plugin if the  =htmlize.el= file is not present.

    local url="https://github.com/hniksic/emacs-htmlize"
    local repo="emacs-htmlize"

    [[ -f ${HTMLIZE} ]] || (
        cd ${DOCS}/html
        ${srcdir}/missing git clone ${url} \
            && cp ${repo}/htmlize.el ${HTMLIZE} \
            && rm -rf ${repo}
        cd -
    )

    # Assert htmlize is installed
    [[ -f ${HTMLIZE} ]] \
        || exit 1
}



function extract_doc()
{
#   Extracts documentation from an org-mode file.

    local org=$1
    local local_html=${ORG}/${org%.org}.html
    local local_text=${ORG}/${org%.org}.txt
    local html=${DOCS}/html/${org%.org}.html
    local text=${DOCS}/text/${org%.org}.txt

    if [[ -f ${html} && ${org} -ot ${html} ]]
    then
        return
    fi
    ${srcdir}/missing emacs --batch  \
          --load ${HTMLIZE}          \
          --load ${CONFIG_DOC}       \
          ${org}                     \
          --load ${CONFIG_TANGLE}    \
          -f org-html-export-to-html \
          -f org-ascii-export-to-ascii
    mv ${local_html} ${DOCS}/html
    mv ${local_text} ${DOCS}/text

}




function main() {

    # Install htmlize if needed
    install_htmlize || exit 2

    # Create documentation
    cd ${ORG} \
        || exit 3

    for i in *.org
    do
        echo "=======  ${i} ======="
        extract_doc ${i}
    done

    if [[ $? -eq 0 ]]
    then
        cd ${DOCS}/html
        rm -f index.html
        ln README.html index.html
        exit 0
    else
        exit 3
    fi

}
main
