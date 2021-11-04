#!/bin/bash

if [ -d throwaway ]
then
  echo used stuff already installed!
else
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt-get update
    sude apt-get install libboost-regex-dev
    sudo apt-get install libboost-dev libboost-test-dev libboost-filesystem-dev libboost-program-options-dev libboost-graph-parallel-dev python g++ libtbb-dev libeigen3-dev python-simplejson libgoogle-perftools-dev install libboost-graph-dev
    #sudo apt-get install libfabric1:i386 libpsm-infinipath1:i386 libhwloc-plugins:i386
    sudo apt-get install mpi openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.10 libopenmpi-dev
    sudo apt-get install libboost-all-dev
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    #brew install boost
    #brew install openmpi
    brew upgrade python
    sudo easy_install simplejson
    sudo easy_install libboost-dev libboost-test-dev libboost-filesystem-dev libboost-program-options-dev libboost-graph-parallel-dev python g++ libtbb-dev libeigen3-dev python-simplejson libgoogle-perftools-dev
    sudo easy_install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.10 libopenmpi-dev

  fi

  mkdir throwaway
fi

if [ -d sferes2 ]
then
  echo nothing to pull
else
  git clone https://github.com/sferes2/sferes2.git
fi


cd sferes2

pwd

./waf configure
./waf configure --mpi /usr/lib/openmpi-bin/
./waf configure --mpi /usr/lib/x86_64-linux-gnu/openmpi
#Edit the ea.hpp file
cp ../ea_edited.hpp sferes/ea/ea.hpp
#MAYBE THIS NEED A PREVIOUS ./waf build ?
./waf build
./waf --tests
#Now we create the different executables
#Use different bashes

pwd


printf "initialisation\n"
for N in 1 #{1..5..1}
do
  FILE=func${N}
  ./waf --create $FILE
  cp ../good_wscript exp/$FILE/wscript
  sed -i "s/_FUNCNAME_/$FILE/" exp/$FILE/wscript
  #CHANGE FITNESS FUNCTION HERE
  cp ../$FILE.cpp exp/$FILE/$FILE.cpp
  ./waf --exp $FILE
done

printf "Computation start\n"
#running:
RESULT_F=results.dat
for N in 1 #{1..5..1}
do
  FILE=func${N}
  cd build/exp/$FILE
  rm $RESULT_F
  touch $RESULT_F
  for T in {1..16..1} #number of threads
  do
    printf "Numberofthreads: ${T}\n"
    for iterations in {1..10..1} # as often as we do PAPSO
    do
      printf "Numberofcalls: ${iterations}\n"
      ./$FILE  #HERE edit number of cores
      #timing in timing.dat
      printf "${N} ${T} " >> $RESULT_F
      cat timing.dat >> $RESULT_F
      printf " " >> $RESULT_F
      #best values in subfolder with weird name in bestfit.dat
      cd timingfolder
      echo $(tail -n 1 bestfit.dat)  | rev | cut -d' ' -f 1 | rev >> ../$RESULT_F
      cd ../
      #printf "\n" >> $RESULT_F
    done
  done
  cd ../../../
done

printf "\nFinished\n"
