#!/bin/bash

OUTPUT=$1

for i in README.org \
             qmckl.org \
             qmckl_context.org \
             qmckl_error.org \
             qmckl_precision.org \
             qmckl_memory.org \
             qmckl_distance.org \
             qmckl_ao.org \
             qmckl_footer.org \
             test_qmckl.org
do
    cat $i >> $1
done
