#!/bin/bash

#
# convert allele codes from affycodes into numerical codes
# affycodes assumed to be phased therefore AB =/= BA
#
#

inpfile=$1

head -n 1 ${inpfile}

#convert hkxhk into separate maternal and paternal markers
tail -q -n +2 ${inpfile} | grep '<hkxhk>' | sed "s/,A./,0/g; s/,B./,1/g; s/Affx-/MAffx-/g"
tail -q -n +2 ${inpfile} | grep '<hkxhk>' | sed "s/,.A/,0/g; s/,.B/,1/g; s/Affx-/PAffx-/g"

#convert lmxll and nnxnp
tail -q -n +2 ${inpfile} | grep -v '<hkxhk>' | sed "s/,AA/,0/g; s/,AB/,1/g; s/,BA/,1/g; s/,BB/,2/g"
