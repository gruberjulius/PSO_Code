#!/bin/bash

i=0
for subdir in DatenSafe/conventionSafe
do
    for file in ${subdir}/*; do
        i=$((i+1))
        touch tmp.txt
        destdir=tmp.txt
        if [ -f "$destdir" ]
        then 
            echo "$file" >> "$destdir"
        fi
        python3 preperationfordata.py < ${destdir}
        rm ${destdir}
    done
    for file in ${subdir}/*/*; do
        i=$((i+1))
        touch tmp.txt
        destdir=tmp.txt
        if [ -f "$destdir" ]
        then 
            echo "$file" >> "$destdir"
        fi
        python3 preperationfordata.py < ${destdir}
        rm ${destdir}
    done
    for file in ${subdir}/*/*/*; do
        i=$((i+1))
        touch tmp.txt
        destdir=tmp.txt
        if [ -f "$destdir" ]
        then 
            echo "$file" >> "$destdir"
        fi
        python3 preperationfordata.py < ${destdir}
        rm ${destdir}
    done
    
done
