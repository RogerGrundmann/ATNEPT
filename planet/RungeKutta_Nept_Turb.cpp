#include "cNeptuneModel.h"
#include "TbNept.h"
#include "CwbNept.h"        // ATNEPT_CWB_DIAG, print-only
#include "UbNept.h"         // ATNEPT_UB_DIAG, print-only
#include "VbNept.h"         // ATNEPT_VBUDGET, print-only

using namespace std;

/*
 * Area-weighted horizontal mean, at each level, of the expression the buoyancy takes the anomaly
 * of. Mirrors ATSAT's computeBuoyancyRefLevel(), which mirrors ATJUP's.
 *
 * The mean has to be the mean OF the quantity whose anomaly is taken, or the anomaly no longer has
 * zero mean at that height — which is the whole point: with it subtracted, only horizontal density
 * contrasts drive vertical motion and hydrostatic balance is left to carry the mean. sin(theta) is
 * the spherical area weight; cells with no positive temperature are skipped, which also catches
 * NaN, and a non-finite contribution is dropped rather than poisoning the level.
 */
void cNeptuneModel::computeBuoyancyRefLevel(){
    if((int)buoy_ref_level.size() != im) buoy_ref_level.assign(im, 0.0);

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double sum = 0.0, wsum = 0.0;
        for(int j = 0; j < jm; j++){
            const double wgt = sin(the.z[j]);              // spherical area weight
            for(int k = 0; k < km; k++){
                if(!(t.x[i][j][k] > 0.0)) continue;        // also catches NaN
                const double b = g * p_stat.x[i][j][k]
                               / (r_mix * R_mix * t.x[i][j][k] * t_ref);
                if(!std::isfinite(b)) continue;
                sum  += wgt * b;
                wsum += wgt;
            }
        }
        buoy_ref_level[i] = (wsum > 0.0) ? sum / wsum : 0.0;
    }
}

/*
 * Hydrostatic pressure perturbation: p_hydro(r) = integral of the buoyancy ANOMALY from the base,
 * so d(p_hydro)/dr = buoyancy by construction. Ported from ATJUP's computeHydrostaticPressure().
 *
 * THE PROBLEM IT SOLVES, measured on this model rather than assumed. The buoyancy term in rhs_u is
 * ~1e5 too small because p_stat is in bar where the ideal-gas density needs pascals. Restoring that
 * factor was measured over 224 iterations and moved the answer 0.5 % — from OLR/in 30.093 to
 * 30.247 — because the pressure projection absorbs a larger radial body force and returns a
 * matching dpdr. The balance is enforced numerically either way, so nothing is left over to drive
 * vertical motion, and this model's radial velocity sits at ~1e-4 where ATSAT's is ~1e-1. With no
 * overturning, thermal diffusion flattens the column unopposed.
 *
 * Splitting it off analytically is what changes that: the radial force balances exactly and by
 * construction, and what enters the momentum equation is the HORIZONTAL gradient of p_hydro —
 * smaller by the 1/r that the horizontal derivative carries, and the part that physically drives a
 * circulation. p_dyn is left with the barotropic and non-hydrostatic remainder, which is the part
 * it can actually represent.
 *
 * The base is the DEEP boundary: the reference isobaric surface, where the gas is densest and
 * horizontal pressure contrasts are hardest to sustain, while the model top is an arbitrary cut
 * through a continuing atmosphere. ATNEPT_HYDRO_REF=1 integrates downward from the top instead so
 * the choice can be measured; it is not the intended configuration.
 */
void cNeptuneModel::computeHydrostaticPressure(){
    static const bool ref_top = [](){
        const char* e = getenv("ATNEPT_HYDRO_REF"); return e && atoi(e) != 0; }();
    // The same factor rhs_u puts on the buoyancy, read the same way, so the two cannot drift apart.
    static const double buoy_scale_local = [](){
        const char* e = getenv("ATNEPT_BUOY_SCALE"); return e ? atof(e) : 1.0; }();

    const double nd = 1.0e5 * L_atm / (u_0 * u_0);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            auto buoy = [&](int i)->double{
                const double T = t.x[i][j][k];
                if(!(T > 0.0)) return 0.0;
                const double b = g * p_stat.x[i][j][k] / (r_mix * R_mix * T * t_ref)
                               - buoy_ref_level[i];
                const double f = -nd * buoy_scale_local * buoyancy * b;
                return std::isfinite(f) ? f : 0.0;
            };
            // ATNEPT_METRIC_SPLIT=1: the column integral carries the metric term (v2 + w2)/r
            // INSTEAD, as rhs_u forms it. The pressure solver cannot hold the gradient that
            // balances a radial body force -- measured here with ATNEPT_UB_DIAG (the balanced
            // initial p_dyn is 0.39 / 2.73 / 0.19 e-3 m/s per iteration at levels 1 / 20 / 39
            // after 224 iterations, against a metric term of 3.51 / 3.96 / 4.22), and measured
            // four ways in ATOM_Precipitation, where neither a line-implicit solve nor a
            // projection in the loop cured it and the column integral did. With the knob rhs_u
            // drops the term (it and -dp_hydro/dr cancel by construction), and rhs_v / rhs_w
            // take the horizontal gradient of p_hydro, which is smaller by the shell's aspect
            // ratio. The buoyancy keeps its own switch (ATNEPT_HYDRO_SPLIT) and is not added
            // here: on this model it is ~1e5 too small to matter either way.
            if(metric_split()){
                // ATNEPT_CORIOLIS_SPLIT: plus the radial Coriolis force as rhs_u forms it,
                // -Coriolis * (L_atm/u_0) * (-2 omega sin(theta) w) -- with the TRUE sine (the
                // model's sinthe_min() is 0 by default), not the metric floor 0.4 the RHS's table
                // carries: a metric floor has no business in a body force, and since the term
                // leaves rhs_u altogether nothing has to match it there.
                const double cor_ = coriolis_split()
                    ? Coriolis * coriolis_scale() * 2.0 * omega
                      * std::max(sinthe_min(), std::abs(sin(the.z[j]))) : 0.0;
                auto met = [&](int i)->double{
                    return (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k])
                         / metricRadius(rad.z[i]) + cor_ * w.x[i][j][k]; };
                p_hydro.x[0][j][k] = 0.0;
                for(int i = 1; i < im; i++)
                    p_hydro.x[i][j][k] = p_hydro.x[i-1][j][k] + 0.5 * (met(i-1) + met(i)) * dr;
                // ATNEPT_PDYN_GEO=2: the frozen baroclinic pressure of the balanced start
                if(!p_baro.empty())
                    for(int i = 0; i < im; i++) p_hydro.x[i][j][k] += p_baro[((size_t)i * jm + j) * km + k];
                continue;
            }
            if(!ref_top){
                p_hydro.x[0][j][k] = 0.0;
                for(int i = 1; i < im; i++)
                    p_hydro.x[i][j][k] = p_hydro.x[i-1][j][k]
                                       + 0.5 * (buoy(i-1) + buoy(i)) * dr;
            } else {
                p_hydro.x[im-1][j][k] = 0.0;
                for(int i = im-2; i >= 0; i--)
                    p_hydro.x[i][j][k] = p_hydro.x[i+1][j][k]
                                       - 0.5 * (buoy(i+1) + buoy(i)) * dr;
            }
        }
    }
}

void cNeptuneModel::RungeKuttaNept(){
    cout << endl << "      ATNEPT: RungeKuttaNept" << endl;


    auto begin = std::chrono::high_resolution_clock::now();

    // The buoyancy base state and the hydrostatic pressure built from it, refreshed once per RK4
    // step before any stage reads them. Both are inert unless ATNEPT_HYDRO_SPLIT is set.
    computeBuoyancyRefLevel();
    computeHydrostaticPressure();

    // Precompute sin/cos tables — depend only on j.
    // sinthe is clamped to a minimum to prevent 1/sin²θ blow-up near the poles.
    static const double sinthe_min = [](){ const char* e = getenv("ATNEPT_SINTHE_FLOOR"); return e ? atof(e) : 0.4; }();   // ATNEPT_SINTHE_FLOOR (2026-10-12): the hardcoded polar metric floor as a knob, default 0.4 = unchanged
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = std::max(sinthe_min, std::abs(sin(the.z[j])));
        costhe_tbl[j] = cos(the.z[j]);
    }

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr   * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    // Turbulence clamps, used only when the closure is on. k* has a physical ceiling; dis* has a
    // FLOOR because it appears in denominators throughout the closure (nue = k/dis among them) and
    // a zero there is an infinity one step later.
    const bool turb_on_rk = turb_active;
    const double tke_max_nd = 1000.0 / (u_0 * u_0);   // 1000 m2/s2
    constexpr double dis_min_nd = 1.0e-10;            // matches the closure's dis_min

    // Physical bounds from the model, nondimensionalised here — see cNeptuneModel.h for why
    // these are no longer bare literals and what the old ceiling was truncating.
    const double t_min = t_min_K() / t_ref;
    const double t_max = t_max_K() / t_ref;

    // ATNEPT_BC_PHI_PERIODIC (ATPhys.h): with 1 the seam meridian k = 0 is integrated like every
    // other cell (RHSNept takes its phi neighbours round the ring); k = km-1 is its copy, written
    // by bcPhi after the step. With 0 the loops run k = 1..km-2, as before.
    const int k_first = ATPhys::phi_first(ATPhys::phi_periodic<cNeptuneModel>());

    /*
     * FOUR BARRIER-SEPARATED STAGES, which is the whole point of this file.
     *
     * What was here ran all four RK4 stages inside ONE parallel loop over cells, holding k1..k4 as
     * per-thread scalars and overwriting t,u,v,w at the cell as it went — while RHSNept
     * differentiates those same live fields at i+-1, j+-1, k+-1, cells other threads are writing
     * at the same moment. Two runs of the same binary at 24 threads differed in 4 of 7 output
     * files; that race is why.
     *
     * Each stage is now two passes with a barrier between them, and the running sum
     * k1 + 2k2 + 2k3 + k4 lives in the acc_* arrays rather than in registers:
     *
     *   pass A   evaluate RHSNept over every cell, reading a state no one is writing
     *   pass B   fold rhs_* into acc_*, then form the NEXT stage's input from the n-copies
     *
     * c_in is the offset of the next stage's input from y_n, wgt this stage's weight in the sum.
     * Stage 3 forms no input, which is why the write is guarded. This is ATSAT's 71082e7 applied
     * to Neptune, and the acceptance test is its: reproducibility at any thread count.
     */
    for(int stage = 0; stage < 4; stage++){

        const double c_in = (stage == 0 || stage == 1) ? 0.5 * dt : (stage == 2 ? dt : 0.0);
        const double wgt  = (stage == 0 || stage == 3) ? 1.0 : 2.0;
        NeptTb::st().stage_wgt = wgt;                 // ATNEPT_TB_DIAG: weight of the tendencies RHSNept is about to sum
        NeptCwb::st().stage_wgt = wgt;                // ATNEPT_CWB_DIAG: the same weight
        NeptUb::st().stage_wgt = wgt;                 // ATNEPT_UB_DIAG: the same weight
        NeptVb::st().stage = stage;                   // ATNEPT_VBUDGET: records stage 0 only

        // ---- pass A: evaluate the right-hand sides ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 1; j < jm-1; j++){
            CellGeometry geo;
            // metricRadius(): identity unless ATNEPT_METRIC_RADIUS is set, so this is inert by
            // default. See the note at metricRadius() in the model header — this model has no
            // metric radius, i.e. rad.z runs 1..2 and the metric puts the planet's surface
            // L_atm from its centre instead of R, making every HORIZONTAL derivative
            // R/L_atm times too large. ATJUP measured the same defect on itself as 499x and
            // fixed it by shifting rad.z at initialisation; ATSAT routes it through this
            // accessor, as here.
            geo.rm      = metricRadius(rad.z[i]);
            geo.rm2     = geo.rm * geo.rm;
            geo.sinthe  = sinthe_tbl[j];
            geo.sinthe2 = geo.sinthe * geo.sinthe;
            geo.costhe  = costhe_tbl[j];
            geo.cotanthe             = geo.costhe / geo.sinthe;
            geo.inv_rm               = 1.0 / geo.rm;
            geo.inv_rm2              = 1.0 / geo.rm2;
            geo.inv_rmsinthe         = 1.0 / (geo.rm * geo.sinthe);
            geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
            geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
            geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
            geo.inv_2dr   = inv_2dr;
            geo.inv_2dthe = inv_2dthe;
            geo.inv_2dphi = inv_2dphi;
            geo.inv_dr2   = inv_dr2;
            geo.inv_dthe2 = inv_dthe2;
            geo.inv_dphi2 = inv_dphi2;

                for(int k = k_first; k < km-1; k++)
                    cNeptuneModel::RHSNept(i, j, k, geo);
            }
        }

        // ---- pass B: fold into the accumulator, then form the next stage's input ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 1; j < jm-1; j++){
                for(int k = k_first; k < km-1; k++){
                    if(turb_on_rk){
                        acc_tke.x[i][j][k] = (stage == 0) ? wgt * rhs_tke.x[i][j][k]
                                           : acc_tke.x[i][j][k] + wgt * rhs_tke.x[i][j][k];
                        acc_dis.x[i][j][k] = (stage == 0) ? wgt * rhs_dis.x[i][j][k]
                                           : acc_dis.x[i][j][k] + wgt * rhs_dis.x[i][j][k];
                    }
                    acc_t.x[i][j][k] = (stage == 0) ? wgt * rhs_t.x[i][j][k]
                                       : acc_t.x[i][j][k] + wgt * rhs_t.x[i][j][k];
                    acc_u.x[i][j][k] = (stage == 0) ? wgt * rhs_u.x[i][j][k]
                                       : acc_u.x[i][j][k] + wgt * rhs_u.x[i][j][k];
                    acc_v.x[i][j][k] = (stage == 0) ? wgt * rhs_v.x[i][j][k]
                                       : acc_v.x[i][j][k] + wgt * rhs_v.x[i][j][k];
                    acc_w.x[i][j][k] = (stage == 0) ? wgt * rhs_w.x[i][j][k]
                                       : acc_w.x[i][j][k] + wgt * rhs_w.x[i][j][k];
                    acc_ch4.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4.x[i][j][k]
                                       : acc_ch4.x[i][j][k] + wgt * rhs_ch4.x[i][j][k];
                    acc_ch4_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_cloud.x[i][j][k]
                                       : acc_ch4_cloud.x[i][j][k] + wgt * rhs_ch4_cloud.x[i][j][k];
                    acc_ch4_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_ice.x[i][j][k]
                                       : acc_ch4_ice.x[i][j][k] + wgt * rhs_ch4_ice.x[i][j][k];
                    acc_h2o.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o.x[i][j][k]
                                       : acc_h2o.x[i][j][k] + wgt * rhs_h2o.x[i][j][k];
                    acc_h2o_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_cloud.x[i][j][k]
                                       : acc_h2o_cloud.x[i][j][k] + wgt * rhs_h2o_cloud.x[i][j][k];
                    acc_h2o_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_ice.x[i][j][k]
                                       : acc_h2o_ice.x[i][j][k] + wgt * rhs_h2o_ice.x[i][j][k];
                    acc_h2s.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s.x[i][j][k]
                                       : acc_h2s.x[i][j][k] + wgt * rhs_h2s.x[i][j][k];
                    acc_h2s_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s_cloud.x[i][j][k]
                                       : acc_h2s_cloud.x[i][j][k] + wgt * rhs_h2s_cloud.x[i][j][k];
                    acc_h2s_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s_ice.x[i][j][k]
                                       : acc_h2s_ice.x[i][j][k] + wgt * rhs_h2s_ice.x[i][j][k];
                    acc_nh3.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3.x[i][j][k]
                                       : acc_nh3.x[i][j][k] + wgt * rhs_nh3.x[i][j][k];
                    acc_nh3_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_cloud.x[i][j][k]
                                       : acc_nh3_cloud.x[i][j][k] + wgt * rhs_nh3_cloud.x[i][j][k];
                    acc_nh3_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_ice.x[i][j][k]
                                       : acc_nh3_ice.x[i][j][k] + wgt * rhs_nh3_ice.x[i][j][k];
                    acc_nh4sh.x[i][j][k] = (stage == 0) ? wgt * rhs_nh4sh.x[i][j][k]
                                       : acc_nh4sh.x[i][j][k] + wgt * rhs_nh4sh.x[i][j][k];

                    if(stage < 3){
                        t.x[i][j][k] = tn.x[i][j][k] + c_in * rhs_t.x[i][j][k];
                        if(turb_on_rk){
                            tke.x[i][j][k] = AtomUtils::clamp(tken.x[i][j][k]
                                + c_in * rhs_tke.x[i][j][k], 0.0, tke_max_nd);
                            dis.x[i][j][k] = std::max(disn.x[i][j][k]
                                + c_in * rhs_dis.x[i][j][k], dis_min_nd);
                        }
                        u.x[i][j][k] = un.x[i][j][k] + c_in * rhs_u.x[i][j][k];
                        v.x[i][j][k] = vn.x[i][j][k] + c_in * rhs_v.x[i][j][k];
                        w.x[i][j][k] = wn.x[i][j][k] + c_in * rhs_w.x[i][j][k];
                        ch4.x[i][j][k] = ch4n.x[i][j][k] + c_in * rhs_ch4.x[i][j][k];
                        ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + c_in * rhs_ch4_cloud.x[i][j][k];
                        ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + c_in * rhs_ch4_ice.x[i][j][k];
                        h2o.x[i][j][k] = h2on.x[i][j][k] + c_in * rhs_h2o.x[i][j][k];
                        h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + c_in * rhs_h2o_cloud.x[i][j][k];
                        h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + c_in * rhs_h2o_ice.x[i][j][k];
                        h2s.x[i][j][k] = h2sn.x[i][j][k] + c_in * rhs_h2s.x[i][j][k];
                        h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + c_in * rhs_h2s_cloud.x[i][j][k];
                        h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + c_in * rhs_h2s_ice.x[i][j][k];
                        nh3.x[i][j][k] = nh3n.x[i][j][k] + c_in * rhs_nh3.x[i][j][k];
                        nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + c_in * rhs_nh3_cloud.x[i][j][k];
                        nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + c_in * rhs_nh3_ice.x[i][j][k];
                        nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + c_in * rhs_nh4sh.x[i][j][k];
                    }
                }
            }
        }
    }

    // ---- the step itself: y_{n+1} = y_n + dt/6 * (k1 + 2k2 + 2k3 + k4) ----
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            for(int k = k_first; k < km-1; k++){
                t.x[i][j][k] = AtomUtils::clamp(
                    tn.x[i][j][k] + dt * acc_t.x[i][j][k] / 6.0, t_min, t_max);
                if(turb_on_rk){
                    tke.x[i][j][k] = AtomUtils::clamp(std::max(tken.x[i][j][k]
                        + dt * acc_tke.x[i][j][k] / 6.0, 0.0), 0.0, tke_max_nd);
                    dis.x[i][j][k] = std::max(disn.x[i][j][k]
                        + dt * acc_dis.x[i][j][k] / 6.0, dis_min_nd);
                }
                u.x[i][j][k] = un.x[i][j][k] + dt * acc_u.x[i][j][k] / 6.0;
                v.x[i][j][k] = vn.x[i][j][k] + dt * acc_v.x[i][j][k] / 6.0;
                w.x[i][j][k] = wn.x[i][j][k] + dt * acc_w.x[i][j][k] / 6.0;
                ch4.x[i][j][k] = ch4n.x[i][j][k] + dt * acc_ch4.x[i][j][k] / 6.0;
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + dt * acc_ch4_cloud.x[i][j][k] / 6.0;
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + dt * acc_ch4_ice.x[i][j][k] / 6.0;
                h2o.x[i][j][k] = h2on.x[i][j][k] + dt * acc_h2o.x[i][j][k] / 6.0;
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + dt * acc_h2o_cloud.x[i][j][k] / 6.0;
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + dt * acc_h2o_ice.x[i][j][k] / 6.0;
                h2s.x[i][j][k] = h2sn.x[i][j][k] + dt * acc_h2s.x[i][j][k] / 6.0;
                h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + dt * acc_h2s_cloud.x[i][j][k] / 6.0;
                h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + dt * acc_h2s_ice.x[i][j][k] / 6.0;
                nh3.x[i][j][k] = nh3n.x[i][j][k] + dt * acc_nh3.x[i][j][k] / 6.0;
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + dt * acc_nh3_cloud.x[i][j][k] / 6.0;
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + dt * acc_nh3_ice.x[i][j][k] / 6.0;
                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + dt * acc_nh4sh.x[i][j][k] / 6.0;
            }
        }
    }
}
