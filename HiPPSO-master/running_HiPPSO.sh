#!/bin/bash
#requirements g++ and GMPLIB (just do "brew install gmp")

#necessary libs:
sudo apt-get install make g++ libgmp-dev
#get that dope ass repo
git clone https://github.com/alexander-rass/HiPPSO.git
#now make
cd HiPPSO
make

#executable file is now under /src/high_precision_pso

#youre done

#more will come when i understand the function implementation
