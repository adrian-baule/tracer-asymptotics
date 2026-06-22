#include "simulation.h"
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <cstdlib>

SimParams parse_args(int argc, char* argv[]) {
    SimParams p;
    // Defaults
    p.aoup     = {1.0, 1.0, 1e-3};
    p.force    = {1.0, 0.5};
    p.sampling = {0.5, 50.0, -2.0};
    p.T        = 1e3;
    p.N_traj   = 100000;
    p.seed     = 42;
    p.output   = "data/output.csv";

    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--tau_c"))   p.aoup.tau_c       = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_A"))     p.aoup.D_A         = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt"))      p.aoup.dt          = atof(argv[++i]);
        else if (!strcmp(argv[i], "--p"))       p.force.p          = atof(argv[++i]);
        else if (!strcmp(argv[i], "--b_min"))   p.force.b_min      = atof(argv[++i]),
                                                p.sampling.b_min   = atof(argv[i]);
        else if (!strcmp(argv[i], "--b_max"))   p.sampling.b_max   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--gamma"))   p.sampling.gamma   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T"))       p.T                = atof(argv[++i]);
        else if (!strcmp(argv[i], "--N_traj"))  p.N_traj           = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed"))    p.seed             = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--output"))  p.output           = argv[++i];
        else { std::cerr << "Unknown argument: " << argv[i] << "\n"; exit(1); }
    }
    return p;
}

int main(int argc, char* argv[]) {
    SimParams params = parse_args(argc, argv);
    std::cout << "Running with tau_c=" << params.aoup.tau_c
              << " D_A=" << params.aoup.D_A
              << " p=" << params.force.p
              << " b_min=" << params.force.b_min
              << " b_max=" << params.sampling.b_max
              << " T=" << params.T
              << " N_traj=" << params.N_traj << "\n";
    run_simulation(params);
    return 0;
}
