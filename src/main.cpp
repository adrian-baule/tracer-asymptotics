#include "simulation.h"
#include "swimmer_factory.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <filesystem>

int main(int argc, char* argv[]) {
    SimParams p;
    // Defaults
    p.process  = {};
    p.sampling = {0.5, 50.0, -2.0};
    p.T        = 1e3;
    p.N_traj   = 100000;
    p.seed     = 42;
    p.mu       = 1.0;
    p.output   = "";  // derived from output_dir after arg parsing

#if defined(PROCESS_BMSHORT) || defined(PROCESS_BMSHORT_EXACT)
    p.force = {ForceType::GAUSSIAN, 1.0, 0.01, 1.0, 1.0};
    p.sampling.b_min = 0.01;
    p.sampling.b_max = -1.0;  // sentinel: auto-compute as sqrt(4 * D_bm * T)
#elif defined(PROCESS_BMLONG)
    p.force = {ForceType::COULOMB, 1.0, 0.5, 1.0, 1.0};
    p.sampling.b_min = 0.5;
    p.sampling.b_max = -1.0;  // sentinel: auto-compute as sqrt(4 * D_bm * T)
#elif defined(PROCESS_RTP_EXACT)
    // b_min = 1.0 (not 0.5): caps the dipole force at p/b_min^2 so dt = 1e-3
    // resolves the near-core dynamics, and raises the linear-response scale
    // mu* = v_A*b_min^2/p to 1.0 so all four mu values stay weakly coupled.
    p.force = {ForceType::DIPOLE, 1.0, 1.0, 1.0, 1.0};
    p.sampling.b_min = 1.0;
    p.sampling.b_max = -1.0;  // sentinel: auto-compute from RTP effective diffusion
    p.T              = 1e4;
#else
    p.force = {ForceType::DIPOLE, 1.0, 0.5, 1.0, 1.0};
#endif

    std::string output_dir = "data";  // default output directory

    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--tau_c"))      p.process.tau_c  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_A"))        p.process.D_A    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--v_A"))        p.process.v_A    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_r"))        p.process.D_r    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--omega"))      p.process.omega  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--beta"))       p.process.beta   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tau_0"))      p.process.tau_0  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--D_bm"))       p.process.D_bm   = atof(argv[++i]);
        else if (!strcmp(argv[i], "--V_0"))        p.process.V_0    = atof(argv[++i]),
                                                   p.force.V_0      = atof(argv[i]);
        else if (!strcmp(argv[i], "--sigma"))      p.process.sigma  = atof(argv[++i]),
                                                   p.force.sigma    = atof(argv[i]);
        else if (!strcmp(argv[i], "--mu"))         p.mu             = atof(argv[++i]);
        else if (!strcmp(argv[i], "--mu_list")) {
            p.mu_list.clear();
            std::string spec = argv[++i];
            size_t pos = 0;
            while (pos <= spec.size()) {
                size_t comma = spec.find(',', pos);
                if (comma == std::string::npos) comma = spec.size();
                std::string tok = spec.substr(pos, comma - pos);
                if (!tok.empty()) p.mu_list.push_back(atof(tok.c_str()));
                pos = comma + 1;
            }
        }
        else if (!strcmp(argv[i], "--fd_step"))    p.fd_step        = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt"))         p.process.dt     = atof(argv[++i]);
        else if (!strcmp(argv[i], "--p"))          p.force.p        = atof(argv[++i]);
        else if (!strcmp(argv[i], "--b_min"))      p.force.b_min    = atof(argv[++i]),
                                                   p.sampling.b_min = atof(argv[i]);
        else if (!strcmp(argv[i], "--b_max"))      p.sampling.b_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "--gamma"))      p.sampling.gamma = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T"))          p.T              = atof(argv[++i]);
        else if (!strcmp(argv[i], "--N_traj"))     p.N_traj         = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed"))       p.seed           = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--output_dir")) output_dir       = argv[++i];
        else if (!strcmp(argv[i], "--output"))     p.output         = argv[++i];
        else { std::cerr << "Unknown argument: " << argv[i] << "\n"; return 1; }
    }

    // Ensure output_dir has a trailing slash
    if (!output_dir.empty() && output_dir.back() != '/')
        output_dir += '/';

    // Create the output directory if it doesn't exist
    if (!output_dir.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(output_dir, ec);
        if (ec) {
            std::cerr << "Error: cannot create output directory '" << output_dir
                      << "': " << ec.message() << "\n";
            return 1;
        }
    }

    // Derive output path from directory unless --output was given explicitly
    if (p.output.empty())
        p.output = output_dir + "samples_" + process_name() + ".csv";

#if defined(PROCESS_BMSHORT) || defined(PROCESS_BMLONG) || defined(PROCESS_BMSHORT_EXACT)
    if (p.sampling.b_max < 0.0)
        p.sampling.b_max = 2.0 * std::sqrt(4.0 * p.process.D_bm * p.T);
#elif defined(PROCESS_RTP_EXACT)
    // 2D RTP effective diffusion: D_eff = v_A^2 / (d * omega) with d = 2
    if (p.sampling.b_max < 0.0) {
        double D_eff = p.process.v_A * p.process.v_A / (2.0 * p.process.omega);
        p.sampling.b_max = 2.0 * std::sqrt(4.0 * D_eff * p.T);
    }
#endif

#if defined(PROCESS_RTP_EXACT) || defined(PROCESS_BMSHORT_EXACT)
    if (p.mu_list.empty()) {
  #if defined(PROCESS_RTP_EXACT)
        p.mu_list = {0.1, 0.3, 1.0, 3.0};
  #else
        p.mu_list = {0.1, 0.3, 1.0, 3.0};
  #endif
    }
    if (static_cast<int>(p.mu_list.size()) > MAX_MU) {
        std::cerr << "Error: --mu_list has " << p.mu_list.size()
                  << " entries, at most " << MAX_MU << " are supported.\n";
        return 1;
    }
    for (double m : p.mu_list) {
        if (!(m > 0.0)) {
            std::cerr << "Error: --mu_list entries must be positive (got " << m << ").\n";
            return 1;
        }
    }
#endif

    if (process_name() == "levy2" && p.process.tau_0 < 100.0 * p.process.dt)
        std::cerr << "Warning: tau_0/dt = " << p.process.tau_0 / p.process.dt
                  << ". Recommend tau_0 >> dt for levy2.\n";

    std::cout << "Process: " << process_name() << "\n";

#if defined(PROCESS_BMSHORT) || defined(PROCESS_BMSHORT_EXACT)
    std::cout << "D_bm=" << p.process.D_bm << " V_0=" << p.process.V_0
              << " sigma=" << p.process.sigma << "\n";
    std::cout << "b_max = " << p.sampling.b_max << " = 2*sqrt(4*D_bm*T)\n";
#elif defined(PROCESS_BMLONG)
    std::cout << "D_bm=" << p.process.D_bm << " sigma=" << p.force.sigma
              << " b_min=" << p.force.b_min << " mu=" << p.mu << "\n";
    std::cout << "Warning: b_max = " << p.sampling.b_max
              << " = 2*sqrt(4*D_bm*T); scales as sqrt(T).\n";
#elif defined(PROCESS_RTP_EXACT)
    std::cout << "v_A=" << p.process.v_A << " omega=" << p.process.omega
              << " p=" << p.force.p << " b_min=" << p.force.b_min << "\n";
    std::cout << "D_eff = v_A^2/(2*omega) = "
              << p.process.v_A * p.process.v_A / (2.0 * p.process.omega)
              << ";  b_max = " << p.sampling.b_max << " = 2*sqrt(4*D_eff*T)\n";
#else
    std::cout << "tau_c=" << p.process.tau_c << " D_A=" << p.process.D_A
              << " v_A=" << p.process.v_A   << " D_r=" << p.process.D_r
              << " omega=" << p.process.omega << " beta=" << p.process.beta << "\n";
#endif
    std::cout << "dt=" << p.process.dt
              << " b_min=" << p.force.b_min << " b_max=" << p.sampling.b_max
              << " T=" << p.T << " N_traj=" << p.N_traj << "\n";
    std::cout << "Output dir: " << output_dir << "\n";

    run_simulation(p);
    return 0;
}
