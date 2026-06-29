#pragma once
#include <string>
#include "force.h"
#include "sampling.h"
#include "swimmer_factory.h"

struct SimParams {
    ProcessParams process;
    ForceParams force;
    SamplingParams sampling;
    double T;           // total simulation time per trajectory
    int N_traj;         // number of trajectories
    int seed;           // base RNG seed
    double mu = 1.0;    // tracer mobility (used by bmlong for A1/A2 scaling)
    std::string output; // output CSV filename
};

// Run the full simulation, write results to output file
// Each row of output: b, A, w (impact parameter, scattering increment, importance weight)
void run_simulation(const SimParams& params);
