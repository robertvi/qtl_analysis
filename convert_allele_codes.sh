#!/bin/bash

#
# convert allele codes from affycodes (AA/AB/BA/BB) into numerical codes
# AA->0 AB/BA->1 BB->2

inpfile=$1

head -n 1 ${inpfile}
tail -q -n +2 ${inpfile} | sed "s/,AA/,0/g; s/,AB/,1/g; s/,BA/,1/g; s/,BB/,2/g"
