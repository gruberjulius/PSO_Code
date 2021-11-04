#!/bin/bash
#mdatom executive in folder
#coords.inp and params.inp in folder

DIR=outputs
rm -r $DIR
mkdir $DIR

cd $DIR
N=100
for T in 275 #bzw: {init..end..incr}
do
  for S in {10001..200000..100}
  do
    PARAM_F=param${N}_${S}_${T}.inp
    sed "s/_NumberAtoms_/${N}/; s/_TargetTemperature_/${T}/; s/_MDSteps_/${S}/" ../params_S.inp > $PARAM_F
    OUTPUT_F=result.out
    ./../mdatom $PARAM_F > $OUTPUT_F

  done
done
python3 ../plot_end2end_S.py

cd ../
