#pragma once
#include <string>
#include "aoup.h"
#include "force.h"
#include "sampling.h"

struct SimParams {
    AOUPParams aoup;
    ForceParams force;
    SamplingParams sampling;
    double T;           // total simulation time per trajectory
    int N_traj;         // number of trajectories
    int seed;           // base RNG seed
    std::string output; // output CSV filename
};

// Parse command line arguments into SimParams
SimParams parse_args(int argc, char* argv[]);

// Run the full simulation, write results to output file
// Each row of output: b, A, w (impact parameter, scattering increment, importance weight)
void run_simulation(const SimParams& params);
