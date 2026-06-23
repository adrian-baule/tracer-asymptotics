#include "simulation.h"
#include "swimmer_factory.h"
#include <iostream>
#include <cstring>
#include <cstdlib>

int main(int argc, char* argv[]) {
    SimParams p;
    // Defaults
    p.process  = {};                       // ProcessParams defaults
    p.force    = {1.0, 0.5};
    p.sampling = {0.5, 50.0, -2.0};
    p.T        = 1e3;
    p.N_traj   = 100000;
    p.seed     = 42;
    p.output   = "data/output_" + process_name() + ".csv";

    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--tau_c"))  p.process.tau_c  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_A"))    p.process.D_A    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--v_A"))    p.process.v_A    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_r"))    p.process.D_r    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--omega"))  p.process.omega  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--beta"))   p.process.beta   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tau_0"))  p.process.tau_0  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt"))     p.process.dt     = atof(argv[++i]);
        else if (!strcmp(argv[i], "--p"))      p.force.p        = atof(argv[++i]);
        else if (!strcmp(argv[i], "--b_min"))  p.force.b_min    = atof(argv[++i]),
                                               p.sampling.b_min = atof(argv[i]);
        else if (!strcmp(argv[i], "--b_max"))  p.sampling.b_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "--gamma"))  p.sampling.gamma = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T"))      p.T              = atof(argv[++i]);
        else if (!strcmp(argv[i], "--N_traj")) p.N_traj         = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed"))   p.seed           = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--output")) p.output         = argv[++i];
        else { std::cerr << "Unknown argument: " << argv[i] << "\n"; return 1; }
    }

    if (process_name() == "levy2" && p.process.tau_0 < 100.0 * p.process.dt)
        std::cerr << "Warning: tau_0/dt = " << p.process.tau_0 / p.process.dt
                  << ". Recommend tau_0 >> dt for levy2.\n";

    std::cout << "Process: " << process_name() << "\n"
              << "tau_c=" << p.process.tau_c << " D_A=" << p.process.D_A
              << " v_A=" << p.process.v_A   << " D_r=" << p.process.D_r
              << " omega=" << p.process.omega << " beta=" << p.process.beta
              << " dt=" << p.process.dt << "\n"
              << "p=" << p.force.p << " b_min=" << p.force.b_min
              << " b_max=" << p.sampling.b_max
              << " T=" << p.T << " N_traj=" << p.N_traj << "\n";

    run_simulation(p);
    return 0;
}
