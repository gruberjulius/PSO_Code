#pragma once

#ifdef _MSC_VER
    #pragma warning(disable : 4267 4244 4100)

    #include <algorithm> // triggers errors for std::min if not included with MSVC (19.0)
#endif

#include "base/util.hpp"
//removed GA for sake of euler (has some compilation error)
//#include "base/algorithms/GeneticAlgorithm.hpp" // Moe.hpp, Mutations.hpp, Crossovers.hpp
#include "base/algorithms/DifferentialEvolution.hpp"
#include "base/algorithms/ParticleSwarm.hpp"
#include "base/algorithms/SimulatedAnnealing.hpp"

//EDIT START
#include "base/algorithms/ParallelSynchronusPS0_kristof.hpp"
#include "base/algorithms/ParallelSynchronusPS0_k_test.hpp"
#include "base/algorithms/ParallelSynchronousPSO_k_Reduce_Bcast.hpp"
#include "base/algorithms/ParallelSynchronousPSO_k_AllReduce.hpp"
#include "base/algorithms/ParallelSynchronousPSO_k_AR_smaller.hpp"
#include "base/algorithms/ParallelAsychronousPSO.hpp"
#include "base/algorithms/ParallelSynchronusPS0.hpp"
#include "base/algorithms/ParallelAsychronousPSO_precompute_vel.hpp"
#include "base/algorithms/ParallelAsynchronousPSO_precompute_vel_nonblocking.hpp"
#include "base/algorithms/ParallelAsychronousPSO_multi_obj.hpp"
#include "base/algorithms/ParallelAsychronousPSO_send_move.hpp"
#include "base/algorithms/ParallelAsychronousPSO_send_move_nnblc.hpp"

#include "base/algorithms/ParallelAsychronousPSO_send_move_smaller_packets.hpp"
#include "base/algorithms/ParallelAsychronousPSO_send_move_part_unroll.hpp"
//#include "base/algorithms/ParallelAsychronousPSO_precompute_vel_2.hpp"
//#include "base/algorithms/evaluatefunction.hpp"
/*
#include "base/algorithms/ParallelAsychronousPSO_send_move_termcrit.hpp"*/
//EDIT END
