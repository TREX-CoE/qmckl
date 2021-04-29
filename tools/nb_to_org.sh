#!/bin/bash
# $ nb_to_org.sh notebook.ipynb
# produces the org-mode file notebook.org

set -e

nb=$(basename $1 .ipynb)
jupyter nbconvert --to markdown ${nb}.ipynb --output ${nb}.md
pandoc ${nb}.md -o ${nb}.org
rm ${nb}.md

