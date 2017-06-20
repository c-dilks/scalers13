#!/bin/bash

touch fitdata.dat
rm fitdata.dat
touch fitdata.dat

for x in {0..40}; do
  root -b -q fit.C'('$x')'
done

cat fitdata.dat
