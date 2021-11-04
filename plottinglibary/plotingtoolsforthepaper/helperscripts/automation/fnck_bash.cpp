#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <vector>


int main(int argc, char * argv[])
{
    // create new scripts for all types of plots
    // create scripts only for specific type of plot:
    // compare, pspso, papso, compiler
    // Possible inputs:
    // all, vs, ps, pa, comp
    std::string what_to_plot = argv[1];
    // time, fitness
    // Possible inputs:
    // c, f
    std::string tf = argv[2];
    // With or without deviation
    std::string d = argv[3];
    int dev = std::atoi(argv[3]);
    // regular, compiler, particle
    // Possible inputs:
    // r, c, p
    std::string rcp = argv[4];

    std::string OUTPUT_FILE = "helperscript";
    std::string BASE = "../DatenSafe/conventionSafe";
    std::vector<std::string> BASE_V;
    std::string TIMES_PSPSO = "/times_pspso";
    std::string TIMES_PAPSO = "/times_papso";
    std::string TIMES_PSO = "/times_pso";
    std::string THE_CHOSEN_ONE;
    std::string PARTICLE;
    std::string COMPILER;
    std::vector<std::string> builds;

    std::string TAIL = "kristof";

    if(rcp == "r"){
        PARTICLE = "";
        //BASE_V.reserve(1);
        //BASE_V = {BASE};
    }
    if(rcp == "c"){
        THE_CHOSEN_ONE = TIMES_PSPSO + "_" + TAIL + ".dat";
        COMPILER = "/compiler_flags";
        builds = {"build_march_o1", "build_march_o2", "build_march_o3", "build_march_x", "build_x_o1", "build_x_o2", "build_x_o3", "build_x_x"};

        BASE_V.reserve(8);
        for(int i = 0; i < builds.size(); ++i){
            std::string N_BASE = BASE + COMPILER + "/" + builds[i];
            BASE_V.push_back(N_BASE);
        }
    }
    if(rcp == "p"){
        PARTICLE = "_4000_500";
        BASE += "/particle_results";
    }

    if(what_to_plot == "vs" || what_to_plot == "all"){
        
        for(int i = 0; i < 6; ++i){
            std::string w_file_name = OUTPUT_FILE;
            w_file_name += "_compare";
            w_file_name += PARTICLE + "_";
            w_file_name += tf + d + "_";
            std::string i_ = std::to_string(i);
            w_file_name += i_;
            w_file_name += ".txt";

            std::ofstream h_script(w_file_name); 
            h_script << "3\n"
                << BASE << TIMES_PSPSO << "_k_AllReduce" << PARTICLE << ".dat\n"
                << BASE << TIMES_PAPSO << "_send" << PARTICLE << ".dat\n"
                << BASE << TIMES_PSO << PARTICLE << ".dat\n";
            if(dev){
                h_script << "d" << "\n";
            }  
            else{
                h_script << tf << "\n";
            }
            h_script << "PSO vs parallel implementations (function " << i_ << ")\n"
                << "../PlotsSafe/Compare/compare" << PARTICLE << "_" << d << tf << "_" << i_ << "\n"
                << "Threads\n";
            
            if(tf == "c"){
                h_script << "Time\n";
            }
            else if(tf == "f"){
                h_script << "Fitness\n";
            }

            h_script << "PSO vs parallel implementations (function " << i_ << ")\n"
                << i_ << "\n"
                << std::to_string(i+1) << "\n"
                << "PSPSO All Reduce (AR)\n"
                << "n\n"
                << "PAPSO Send\n"
                << "n\n"
                << "PSO\n"
                << "y\n";
            h_script.close();
        }
    }
    if(what_to_plot == "ps" || what_to_plot == "all"){
        
        for(int i = 0; i < 6; ++i){
            std::string w_file_name = OUTPUT_FILE;
            w_file_name += "_pspso";
            w_file_name += PARTICLE + "_";
            w_file_name += tf + d + "_";
            std::string i_ = std::to_string(i);
            w_file_name += i_;
            w_file_name += ".txt";

            std::ofstream h_script(w_file_name); 
            h_script << "5\n"
                << BASE << TIMES_PSPSO << "_k_AR_smaller" << PARTICLE << ".dat\n"
                << BASE << TIMES_PSPSO << "_k_AllReduce" << PARTICLE << ".dat\n"
                << BASE << TIMES_PSPSO << "_k_Reduce_Bcast" << PARTICLE << ".dat\n"
                << BASE << TIMES_PSPSO << "_k_test" << PARTICLE << ".dat\n"
                << BASE << TIMES_PSPSO << "_kristof" << PARTICLE << ".dat\n";
            if(dev){
                h_script << "d" << "\n";
            }  
            else{
                h_script << tf << "\n";
            }
            h_script << "PSPSO comparing all versions (function " << i_ << ")\n"
                << "../PlotsSafe/PSPSO/PSPSO_version_comp" << PARTICLE << "_" << d << tf << "_" << i_ << "\n"
                << "Threads\n";

            //std::cout << BASE << TIMES_PSPSO << "_k_AR_smaller" << PARTICLE << ".dat\n";

            if(tf == "c"){
                h_script << "Time\n";
            }
            else if(tf == "f"){
                h_script << "Fitness\n";
            }

            h_script << "PSPSO comparing all versions (function " << i_ << ")\n"
                << i_ << "\n"
                << std::to_string(i+1) << "\n"
                << "PSPSO AR smaller\n"
                << "n\n"
                << "PSPSO All Reduce (AR)\n"
                << "n\n"
                << "PSPSO Reduce\n"
                << "n\n"
                << "PSPSO test\n"
                << "n\n"
                << "PSPSO kristof\n"
                << "n\n";
            h_script.close();
        }
    }
    if(what_to_plot == "pa" || what_to_plot == "all"){
        
        for(int i = 0; i < 6; ++i){
            std::string w_file_name = OUTPUT_FILE;
            w_file_name += "_papso";
            w_file_name += PARTICLE + "_";
            w_file_name += tf + d + "_";
            std::string i_ = std::to_string(i);
            w_file_name += i_;
            w_file_name += ".txt";

            std::ofstream h_script(w_file_name); 
            h_script << "5\n"
                << BASE << TIMES_PAPSO << "_pre" << PARTICLE << ".dat\n"
                << BASE << TIMES_PAPSO << "_pre_nonblc" << PARTICLE << ".dat\n"
                << BASE << TIMES_PAPSO << "_send" << PARTICLE << ".dat\n"
                << BASE << TIMES_PAPSO << "_send_part_unroll(no_2_3_threads)" << PARTICLE << ".dat\n"
                << BASE << TIMES_PAPSO << "_send_smaller_packets" << PARTICLE << ".dat\n";
            if(dev){
                h_script << "d" << "\n";
            }  
            else{
                h_script << tf << "\n";
            }
            h_script << "PAPSO comparing all versions (function " << i_ << ")\n"
                << "../PlotsSafe/PAPSO/PAPSO_version_comp_" << d << tf << "_" << i_ << "\n"
                << "Threads\n";
            
            if(tf == "c"){
                h_script << "Time\n";
            }
            else if(tf == "f"){
                h_script << "Fitness\n";
            }

            h_script << "PAPSO comparing all versions (function " << i_ << ")\n"
                << i_ << "\n"
                << std::to_string(i+1) << "\n"
                << "PAPSO Precompute Velocity\n"
                << "n\n"
                << "PAPSO Precompute Velocity Non Blocking\n"
                << "n\n"
                << "PAPSO Send\n"
                << "n\n"
                << "PAPSO Send Part Unroll\n"
                << "n\n"
                << "PAPSO Send Smaller Packets\n"
                << "n\n";
            h_script.close();
        }
    }
    if((what_to_plot == "comp" || what_to_plot == "all") && rcp == "c")
    {       
        for(int i = 0; i < 6; ++i)
        {
            std::string w_file_name = OUTPUT_FILE;
            w_file_name += "_compiler";
            w_file_name += PARTICLE + "_";
            w_file_name += tf + d + "_";
            std::string i_ = std::to_string(i);
            w_file_name += i_;
            w_file_name += ".txt";

            std::ofstream h_script(w_file_name); 
            h_script << std::to_string(BASE_V.size()) << "\n";
            for(int k = 0; k < BASE_V.size(); ++k){
                h_script << BASE_V[k] << THE_CHOSEN_ONE << "\n";
            };
            if(dev){
                h_script << "d" << "\n";
            }  
            else{
                h_script << tf << "\n";
            }
            h_script << "Comparing compiler flags for " << THE_CHOSEN_ONE << " (function " << i_ << ")\n"
                << "../PlotsSafe/Compiler_flags/" << TAIL << "_" << d << tf << "_" << i_ << "\n"
                << "Threads\n";
            
            if(tf == "c"){
                h_script << "Time\n";
            }
            else if(tf == "f"){
                h_script << "Fitness\n";
            }

            h_script << "Comparing compiler flags for " << THE_CHOSEN_ONE << " (function " << i_ << ")\n"
                << i_ << "\n"
                << std::to_string(i+1) << "\n";

            for(int k = 0; k < BASE_V.size(); ++k){
                h_script << builds[k] << "\n"
                << "n\n";
            }
            h_script.close();
            
        }
    }
}