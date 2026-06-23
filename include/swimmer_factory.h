#pragma once
#include "swimmer.h"
#include "aoup_swimmer.h"
#include "abp_swimmer.h"
#include "rtp_swimmer.h"
#include "levy_swimmer.h"
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
    // Levy
    double beta   = 1.5;
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
#elif defined(PROCESS_LEVY)
    return std::make_unique<LevySwimmer>(p.v_A, p.beta, p.dt);
#else
#error "No process defined. Compile with -DPROCESS_AOUP, -DPROCESS_ABP, -DPROCESS_RTP, or -DPROCESS_LEVY"
#endif
}

inline std::string process_name() {
#if defined(PROCESS_AOUP)
    return "aoup";
#elif defined(PROCESS_ABP)
    return "abp";
#elif defined(PROCESS_RTP)
    return "rtp";
#elif defined(PROCESS_LEVY)
    return "levy";
#else
    return "unknown";
#endif
}
