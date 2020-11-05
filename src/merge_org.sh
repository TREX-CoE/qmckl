#!/bin/bash

for i in README.org \
             qmckl.org \
             qmckl_memory.org \
             qmckl_context.org \
             qmckl_distance.org \
             qmckl_ao.org \
             test_qmckl.org
do
    cat $i >> merged_qmckl.org
done
