

We will need:  
-one bash script for setting up sferes2 on any system, create and compile the executables with our different benchmark functions  
-the main benchmark bash script should call all those different executables nad then analyse the results



# 1: Setting PGA up

In the terminal: in root I think

`sudo apt-get install libboost-dev libboost-test-dev libboost-filesystem-dev libboost-program-options-dev libboost-graph-parallel-dev python g++ libtbb-dev libeigen3-dev python-simplejson libgoogle-perftools-dev`

*For openmpi:*  
`sudo apt-get install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.10 libopenmpi-dev`

`git clone https://github.com/sferes2/sferes2.git`  
`cd sferes2`  
`./waf configure`  
`./waf configure --mpi /usr/lib/openmpi-bin/ or /usr/lib/x86_64-linux-gnu/openmpi for me sometimes`   
 > EDIT the ea.hpp file st it becomes ea_edited.hpp  
`./waf build`  
`./waf --tests`

*They should all be green. Lol there is one red for me oh well*  
To create the different executables:

`./waf --create new_project`

In sferes2/exp/new_project one can edit the new_project.cpp file.    
Don't forget to change the wscript file too.
Consider writing something like this in the wscript:

``from waflib.Configure import conf  
def options(opt):  
    pass  
def build(bld):  
    bld.program(features = 'cxx',  
            source = 'executable_name.cpp',  
            includes = '. ../../',  
            uselib = 'TBB BOOST EIGEN PTHREAD MPI',  
            use = 'sferes2',  
            target = 'executable_name')``


Now we build the new_project:  

`./waf --exp new_project`

In sferes2/build/exp/new_project there now exists an executable with the name given in the wscript file earlier.  
Important for us is the fitness function.  



# 2: Calling he PGA executables and retreiving informations/results. Can be done solely with the executables.

`./executable_name`  
`cd executable_name_[year]-[month]-[day]_[hour]_[minute]_[seconds]_[PID]`  
`cat bestfit.dat`  *: contains for each generation the value of the best fitness.*  
`./../executable_name --load gen_[inde_of_generation] -o output_file_name -n number`

gen_[generation_index] for a generation | output_file_name contains a description of the best individual of this generation (?) | number: how many individuals you want in the output file | additional -s option to specify the statistics number/index in function of your statistics list  
The output_file indicates the 10 parameters that correspond to the best individual of the selected generation.  

You might analyse the data in bestfit.dat with the ploting function I put in this directory.

#3: Extracting timing and best value:  

In the executable_directory you find the timing in a readable form in:
`timing.dat`  
In the result folder you can use a ploting script to show how we converge to the solution or extract the last value in the last line ie best value of the final generation.  


**We could write for each function we want to test a separate executable.**


>NB: I did not find out how to get the "debug" and "default" builds as writen on the documentation.
