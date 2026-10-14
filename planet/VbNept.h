/*
 * Neptune Atmosphere Circulation Model (ATNEPT)
 * ATNEPT_VBUDGET=<n> -- the meridional momentum budget against latitude, PRINT-ONLY, default off;
 * printed every n iterations.
 *
 * Written 2026-10-14 for the continuity residuum of the balanced start (ATNEPT_GEO_BALANCE): it
 * rises through the run and sits poleward of 75 deg at ~185 km, drifting poleward. Ported from
 * ATURAN's ATURAN_VBUDGET (printMeridionalBudget), with two additions: the transport term is split
 * into its advective and its metric part ((u v - w2 cot) / r), and the five parts of the
 * continuity residuum are printed for the same rows, as Reporting.h forms them.
 *
 * One level (ATNEPT_VBUDGET_I, default 37 = 185 km), one meridian (ATNEPT_VBUDGET_K, default 180),
 * every ATNEPT_VBUDGET_JSTEP-th row (default 10) plus the five rows next to each pole. The terms
 * are those of the first Runge-Kutta stage of the last step, nondimensional, with the sign they
 * enter rhs_v. In RHSNept one (i,k) is one thread per row j, so every slot has one writer.
 *
 * The state lives in a function-local static, as in UbNept.h.
 */

#pragma once

#include <vector>

namespace NeptVb {

    // the terms of rhs_v, with the sign they enter it
    enum { T_DP = 0, T_DPH, T_ADV, T_MET, T_DIFF, T_COR, T_CEN, T_TOT, NTERM };

    struct State {
        bool on;                                     // ATNEPT_VBUDGET > 0
        int  stage;                                  // RK stage being evaluated; recorded at 0 only
        int  i, k;                                   // level and meridian
        std::vector<double> term[NTERM];             // [term][j]
    };

    inline State& st(){
        static State s = State();
        return s;
    }
}
