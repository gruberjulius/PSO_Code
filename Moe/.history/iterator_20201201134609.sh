#!/bin/bash

if [[ -f results/times.txt ]]
then
    rm results/times.txt
    touch results/times.txt
fi

for i in 1000 10000 10000 0 1000000
do
  echo "$i " 
  ./ps_functions 
  <i
done

exit 0