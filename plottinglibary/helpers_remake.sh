#!/bin/bash

# to be called from automation folder
# if moved to different location path has to be changed

cd helperscripts/automation

g++ fnck_bash.cpp

for dev in 0 1
do
    for tfd in c f
    do
        for rcp in r c p
        do
            for filename in *.txt; do
                ./a.out all ${tfd} ${dev} ${rcp}
            done
        done
    done
done
