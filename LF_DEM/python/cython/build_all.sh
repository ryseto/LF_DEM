#!/bin/sh

PYTHON=/opt/local/bin/python

for file in setup.*py
do
    $PYTHON $file build_ext --inplace
done
