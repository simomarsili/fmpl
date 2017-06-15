#!/bin/bash
tmp="$(mktemp)"
echo $tmp
for ffile in *.f90;
do
    echo $ffile;
    tail -n+3 $ffile > $tmp;
    cat HEADER.txt $tmp > $ffile;
done
