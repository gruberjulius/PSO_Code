#!/bin/bash

if [[ -f times.txt ]]
then
    rm times.txt
    touch times.txt
fi

for i in 1000 10000 100000 1000000
do
  echo "$i " 
  ./ps_functions 
  <i
done

exit 0