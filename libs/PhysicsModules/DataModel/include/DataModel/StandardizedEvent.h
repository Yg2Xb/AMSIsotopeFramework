#ifndef STANDARDIZED_EVENT_H
#define STANDARDIZED_EVENT_H

#include <bitset>

// The StandardizedEvent struct: The single, unified data object for the analysis layer.
// It contains both simplified/configured variables and direct copies of common raw branches,
// completely decoupling all downstream modules from the raw TTree format.
struct StandardizedEvent {
    // --- General Event Information ---
    unsigned int run;
    unsigned int event;
    bool is_mc;
    bool is_bad_run;
    std::bitset<32> cut_status;

    // --- Standardized & Configured Physics Quantities ---
    float rigidity;         // Rigidity value from the version specified in config.yaml
    float rigidity_error;   // Corresponding error
    float beta_tof;         // TOF beta from the version specified in config.yaml
    
    // --- Directly Mapped Common Raw Branches ---
    // These variables are direct copies from the amstreea for convenience.
    
    // Beta & Kinematics
    float beta_rich_naf;
    float beta_rich_agl;
    float ek_tof;
    float ek_rich_naf;
    float ek_rich_agl;
    float cutoff;

    // Charge
    float charge_l1_unbiased;
    float charge_l1_normal;
    float charge_l2_normal;
    float charge_inner;
    float charge_tof_upper;
    float charge_tof_lower;
    float charge_rich;

    // RICH Quality
    int   rich_pmt;
    int   rich_hit;
    int   rich_used_hits;
    float rich_npe[3];
    bool  rich_is_naf;
    float rich_beta_consistency;
    
    // Tracker Quality
    int   ntrack;
    int   tk_q_hit[2];
    float tk_q_rms[2];
    float tk_inner_q_rms;

    // MC Truth Information (if applicable)
    float generated_rigidity;
    float generated_ek;
    int   particle_id_L1;
    int   particle_id_L2;
};

#endif // STANDARDIZED_EVENT_H