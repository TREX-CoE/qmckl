#!/bin/bash
# Script to build the documentation
#   :PROPERTIES:
#   :header-args:bash: :tangle build_doc.sh :noweb  yes :shebang #!/bin/bash :comments org
#   :END:

#   First define readonly global variables.


readonly DOCS=${top_srcdir}/share/doc/qmckl/
readonly SRC=${top_srcdir}/src/
readonly HTMLIZE=${DOCS}/html/htmlize.el
readonly CONFIG_DOC=${top_srcdir}/tools/config_doc.el
readonly CONFIG_TANGLE=${top_srcdir}/tools/config_tangle.el



# Check that all the defined global variables correspond to files.


function check_preconditions()
{
    if [[ -z ${top_srcdir} ]]
    then
        print "top_srcdir is not defined"
        exit 1
    fi

    for dir in ${DOCS}/html ${DOCS}/text ${SRC}
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
}



# ~install_htmlize~ installs the htmlize Emacs plugin if the
# =htmlize.el= file is not present.


function install_htmlize()
{
    local url="https://github.com/hniksic/emacs-htmlize"
    local repo="emacs-htmlize"

    [[ -f ${HTMLIZE} ]] || (
        cd ${DOCS}/html
        git clone ${url} \
            && cp ${repo}/htmlize.el ${HTMLIZE} \
            && rm -rf ${repo}
        cd -
    )

    # Assert htmlize is installed
    [[ -f ${HTMLIZE} ]] \
        || exit 1
}



# Extract documentation from an org-mode file.


function extract_doc()
{
    local org=$1
    local local_html=${SRC}/${org%.org}.html
    local local_text=${SRC}/${org%.org}.txt
    local html=${DOCS}/html/${org%.org}.html
    local text=${DOCS}/text/${org%.org}.txt

    if [[ -f ${html} && ${org} -ot ${html} ]]
    then
        return
    fi
    emacs --batch                    \
          --load ${HTMLIZE}          \
          --load ${CONFIG_DOC}       \
          ${org}                     \
          --load ${CONFIG_TANGLE}    \
          -f org-html-export-to-html \
          -f org-ascii-export-to-ascii
    mv ${local_html} ${DOCS}/html
    mv ${local_text} ${DOCS}/text

}



# The main function of the script.


function main() {

    check_preconditions || exit 1

    # Install htmlize if needed
    install_htmlize || exit 2

    # Create documentation
    cd ${SRC} \
        || exit 3

    for i in *.org
    do
        echo
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
