cmake_minimum_required(VERSION 2.8)

project(MOE)

option(BUILD_EXAMPLES "Builds examples" ON)
option(BUILD_TESTS "Builds Tests" OFF)
option(DEBUG "Enables debug symbols" OFF)


include_directories( "include" )

if(BUILD_TESTS OR COVERAGE)
    # Add philsquared/Catch
    include_directories( "3rd/catch/single_include" ) # catch.hpp
endif()

SET(CMAKE_CXX_COMPILER mpicxx)#mpicxx

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "mpicxx")

    # using mpicxx
    MESSAGE(STATUS "Using mpicxx!!!")

    add_compile_options(
                    -std=c++14
                    -Wall
                    -pedantic
                    -O3
                    -march=native
                    )
    if(DEBUG)
        #Debug Mode
        add_compile_options(
                    -g
                    )
    else()
        #Release Mode
        add_compile_options(
                #    -O3
                    )
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

    # using Clang
    MESSAGE(STATUS "Using Clang!!!")

    add_compile_options(
                    -std=c++14
                    -Wall
                    -pedantic
                    #-fopenmp
                    )
    if(DEBUG)
        #Debug Mode
        add_compile_options(
                    -g
                    )
    else()
        #Release Mode
        add_compile_options(
                    -O3
                    )
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    # using GCC
    MESSAGE(STATUS "Using GCC!!!")


    add_compile_options(
                    -std=c++14
                    -Wall
                    -pedantic
                    -fopenmp
                    -O3
                    -march=native
                    )
    if(COVERAGE)
        add_compile_options(
                    -g
                    -O0
                    --coverage
                    )
        set(CMAKE_EXE_LINKER_FLAGS
            --coverage
            )
        set(CMAKE_SHARED_LINKER_FLAGS
            --coverage
            )
    elseif(DEBUG)
        #Debug Mode
        add_compile_options(
                    -g
                    )
    else()
        #Release Mode
        add_compile_options(
                    -O3
                    )
    endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

    # using Visual Studio C++

    add_compile_options(/std=c++14
                    /W4
                    )
    if(DEBUG)
        #Debug Mode
    else()
        #Release Mode
        add_compile_options( /O2t )
    endif()
endif()

if(BUILD_EXAMPLES)

    set(SOURCE
        #"examples/paps_functions.cpp"
        "examples/ps_functions.cpp"
        )

    add_executable(psps_functions_kristof examples/psps_functions_kristof.cpp)
    add_executable(psps_functions_k_test examples/psps_functions_k_test.cpp)
    add_executable(psps_functions_k_Reduce_Bcast examples/psps_functions_k_Reduce_Bcast.cpp)
    add_executable(psps_functions_k_AllReduce examples/psps_functions_k_AllReduce.cpp)
    add_executable(psps_functions_k_AR_smaller examples/psps_functions_k_AR_smaller.cpp)

    #add_executable(ps_functions ${SOURCE})
    add_executable(ps_functions examples/ps_functions.cpp)
    add_executable(psps_functions examples/psps_functions.cpp)
    add_executable(paps_functions examples/paps_functions.cpp)
    add_executable(paps_functions_precompute_vel examples/paps_functions_precompute_vel.cpp)
    add_executable(paps_functions_send_move examples/paps_functions_send_move.cpp)
    add_executable(paps_functions_send_move_nnblc examples/paps_functions_send_move_nnblc.cpp)
    add_executable(paps_functions_precompute_vel_nonblocking examples/paps_functions_precompute_vel_nonblocking.cpp)
    add_executable(paps_functions_send_move_smaller_packets examples/paps_functions_send_move_smaller_packets.cpp)
    add_executable(paps_functions_send_move_part_unroll examples/paps_functions_send_move_part_unroll.cpp)

    #add_executable(paps_functions_precompute_vel_2 examples/paps_functions_precompute_vel_2.cpp)
    #add_executable(paps_functions_multi_obj examples/paps_functions_multi_obj.cpp)

    #set(SOURCE
	    #"examples/sa_functions.cpp"
	#        )
#
#    add_executable(sa_functions ${SOURCE})
endif()

if(BUILD_TESTS OR COVERAGE)
    set(SOURCE
        "tests/tests.cpp"
        )

    add_executable(tests ${SOURCE})
endif()
