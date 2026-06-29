#include "simulation.h"
#include "swimmer_factory.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

int main(int argc, char* argv[]) {
    SimParams p;
    // Defaults
    p.process  = {};
    p.sampling = {0.5, 50.0, -2.0};
    p.T        = 1e3;
    p.N_traj   = 100000;
    p.seed     = 42;
    p.mu       = 1.0;
    p.output   = "data/output_" + process_name() + ".csv";

#if defined(PROCESS_BMSHORT)
    // Gaussian force, no singularity so b_min ~ 0; b_max set from T after arg parsing
    p.force = {ForceType::GAUSSIAN, 1.0, 0.01, 1.0, 1.0};
    p.sampling.b_min = 0.01;
    p.sampling.b_max = -1.0;  // sentinel: auto-compute as sqrt(D_bm * T)
#elif defined(PROCESS_BMLONG)
    // Coulomb force: singular at origin, physical hard core b_min = 0.5
    // b_max = sqrt(4 * D_bm * T) — computed after arg parsing
    p.force = {ForceType::COULOMB, 1.0, 0.5, 1.0, 1.0};  // sigma stored in force.sigma
    p.sampling.b_min = 0.5;
    p.sampling.b_max = -1.0;  // sentinel: auto-compute as sqrt(4 * D_bm * T)
#else
    p.force = {ForceType::DIPOLE, 1.0, 0.5, 1.0, 1.0};
#endif

    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--tau_c"))  p.process.tau_c  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_A"))    p.process.D_A    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--v_A"))    p.process.v_A    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_r"))    p.process.D_r    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--omega"))  p.process.omega  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--beta"))   p.process.beta   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tau_0"))  p.process.tau_0  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_bm"))   p.process.D_bm   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--V_0"))    p.process.V_0    = atof(argv[++i]),
                                               p.force.V_0      = atof(argv[i]);
        else if (!strcmp(argv[i], "--sigma"))  p.process.sigma  = atof(argv[++i]),
                                               p.force.sigma    = atof(argv[i]);
        else if (!strcmp(argv[i], "--mu"))     p.mu             = atof(argv[++i]);
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

#if defined(PROCESS_BMSHORT)
    // Auto-set b_max = sqrt(D_bm * T) if the user did not provide --b_max
    if (p.sampling.b_max < 0.0)
        p.sampling.b_max = std::sqrt(p.process.D_bm * p.T);
#elif defined(PROCESS_BMLONG)
    // Auto-set b_max = sqrt(4 * D_bm * T): diffusion length sets the relevant impact
    // parameter range. Trajectories starting beyond this never reach the force region.
    if (p.sampling.b_max < 0.0)
        p.sampling.b_max = std::sqrt(4.0 * p.process.D_bm * p.T);
#endif

    if (process_name() == "levy2" && p.process.tau_0 < 100.0 * p.process.dt)
        std::cerr << "Warning: tau_0/dt = " << p.process.tau_0 / p.process.dt
                  << ". Recommend tau_0 >> dt for levy2.\n";

    std::cout << "Process: " << process_name() << "\n";

#if defined(PROCESS_BMSHORT)
    std::cout << "D_bm=" << p.process.D_bm << " V_0=" << p.process.V_0
              << " sigma=" << p.process.sigma << "\n";
#elif defined(PROCESS_BMLONG)
    std::cout << "D_bm=" << p.process.D_bm << " sigma=" << p.force.sigma
              << " b_min=" << p.force.b_min << " mu=" << p.mu << "\n";
    std::cout << "Warning: b_max = " << p.sampling.b_max
              << " = sqrt(4*D_bm*T); scales as sqrt(T).\n";
#else
    std::cout << "tau_c=" << p.process.tau_c << " D_A=" << p.process.D_A
              << " v_A=" << p.process.v_A   << " D_r=" << p.process.D_r
              << " omega=" << p.process.omega << " beta=" << p.process.beta << "\n";
#endif
    std::cout << "dt=" << p.process.dt
              << " b_min=" << p.force.b_min << " b_max=" << p.sampling.b_max
              << " T=" << p.T << " N_traj=" << p.N_traj << "\n";

    run_simulation(p);
    return 0;
}
