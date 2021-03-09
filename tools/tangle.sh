#!/bin/bash
# Script to tangle the org-mode files
#   :PROPERTIES:
#   :header-args: :tangle tangle.sh :noweb  yes :shebang #!/bin/bash :comments org
#   :END:


# This file was created by tools/Building.org



# This file needs to be run from the QMCKL =src= directory.

# It tangles all the files in the directory. It uses the
# =config_tangle.el= file, which contains information required to
# compute the current file names using for example ~(eval c)~ to get
# the name of the produced C file.

# The file is not tangled if the last modification date of the org
# file is less recent than one of the tangled files.

function tangle()
{
    if [[ -f ${1%.org}.c && $1 -ot ${1%.org}.c ]]
    then return
    elif [[ -f ${1%.org}.f90 && $1 -ot ${1%.org}.f90 ]]
    then return
    fi
    emacs --batch $1 --load=../tools/config_tangle.el -f org-babel-tangle
}


for i in $@
do
        echo "--- $i ----"
        tangle $i
done
