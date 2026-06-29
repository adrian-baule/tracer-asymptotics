#pragma once
#include "swimmer.h"
#include "aoup_swimmer.h"
#include "abp_swimmer.h"
#include "rtp_swimmer.h"
#include "levy1_swimmer.h"
#include "levy2_swimmer.h"
#include "bmshort_swimmer.h"
#include <memory>
#include <string>

struct ProcessParams {
    // AOUP
    double tau_c = 1.0;
    double D_A   = 1.0;
    // ABP / RTP / Levy shared swim speed
    double v_A   = 1.0;
    // ABP
    double D_r   = 1.0;
    // RTP
    double omega  = 1.0;
    // Levy1 and Levy2
    double beta   = 1.5;   // > 1 for levy1, < 1 for levy2
    // Levy2 only
    double tau_0  = 1.0;   // minimum run time, ignored by levy1
    // bmshort swimmer
    double D_bm   = 1.0;   // BM diffusion coefficient (distinct from D_A)
    // bmshort force
    double V_0    = 1.0;   // Gaussian potential amplitude
    double sigma  = 1.0;   // Gaussian force range
    // shared timestep
    double dt     = 1e-3;
};

inline std::unique_ptr<Swimmer> make_swimmer(const ProcessParams& p) {
#if defined(PROCESS_AOUP)
    return std::make_unique<AOUPSwimmer>(p.tau_c, p.D_A, p.dt);
#elif defined(PROCESS_ABP)
    return std::make_unique<ABPSwimmer>(p.v_A, p.D_r, p.dt);
#elif defined(PROCESS_RTP)
    return std::make_unique<RTPSwimmer>(p.v_A, p.omega, p.dt);
#elif defined(PROCESS_LEVY1)
    return std::make_unique<Levy1Swimmer>(p.v_A, p.beta, p.dt);
#elif defined(PROCESS_LEVY2)
    return std::make_unique<Levy2Swimmer>(p.v_A, p.beta, p.dt, p.tau_0);
#elif defined(PROCESS_BMSHORT)
    return std::make_unique<BMShortSwimmer>(p.D_bm, p.dt);
#elif defined(PROCESS_BMLONG)
    return std::make_unique<BMShortSwimmer>(p.D_bm, p.dt);  // same BM dynamics, Coulomb force
#else
#error "No process defined. Compile with -DPROCESS_AOUP, -DPROCESS_ABP, -DPROCESS_RTP, -DPROCESS_LEVY1, -DPROCESS_LEVY2, -DPROCESS_BMSHORT, or -DPROCESS_BMLONG"
#endif
}

inline std::string process_name() {
#if defined(PROCESS_AOUP)
    return "aoup";
#elif defined(PROCESS_ABP)
    return "abp";
#elif defined(PROCESS_RTP)
    return "rtp";
#elif defined(PROCESS_LEVY1)
    return "levy1";
#elif defined(PROCESS_LEVY2)
    return "levy2";
#elif defined(PROCESS_BMSHORT)
    return "bmshort";
#elif defined(PROCESS_BMLONG)
    return "bmlong";
#else
    return "unknown";
#endif
}
