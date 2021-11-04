#!/bin/bash

TESTSETS=$1
ITERATIONS=$2

for TESTID in $(seq 1 $TESTSETS) ; do
    echo Testseries $TESTID
    REFFOLDER="reference_data$TESTID"
    SPECIALOPTIONS="special_options/special_options$TESTID.conf"
    # create clean folder for reference data
    if [ -d $REFFOLDER ] ; then
        rm -r $REFFOLDER
    fi
    mkdir $REFFOLDER

    numoptions=$(cat $SPECIALOPTIONS | wc -l);
    for line in $(seq 1 $numoptions) ;
    do
        tmpconffile=tmp.conf ;
        cp base_configfiles/base_configfile$TESTID.conf $tmpconffile;
        optionDescription=$(head -n$line $SPECIALOPTIONS | tail -n1);
        echo "$optionDescription" | tr '#' '\n' >> $tmpconffile;
        echo "steps $ITERATIONS" >> $tmpconffile
        appendfile=base_configfiles/base_configfile_append$TESTID.conf
        if [ -f "$appendfile" ] ; then
            cat "$appendfile" >> $tmpconffile
        fi
        folder=run$line ;
        mkdir $REFFOLDER/$folder
        cd $REFFOLDER/$folder
        ../../../../../../bin/high_precision_pso c ../../$tmpconffile > stdout.txt 2> stderr.txt;
        STDOUTLINES=$(cat stdout.txt | wc -l)
        if [ $STDOUTLINES -eq 0 ]; then
            rm stdout.txt;
        fi
        STDERRLINES=$(cat stderr.txt | wc -l)
        if [ $STDERRLINES -eq 0 ]; then
            rm stderr.txt;
        fi
        cd ../..
        rm $tmpconffile
    done
done
