#!/bin/bash

OUTPUT=$1

for i in README.org $(cat $QMCKL_ROOT/src/table_of_contents)
do
    cat $i >> $1
done
