/*
 * Atmosphere General Circulation Modell (ATNEPT)
 * Standalone saturation-adjustment and microphysics class for the Neptune model.
 * Declared as friend of cNeptuneModel so it may access all private members
 * through the stored reference.
 *
 * Reference algorithm:
 *   Tao, W.-K., Simpson, J., and McCumber, M.:
 *   "An Ice-Water Saturation Adjustment", AMS Notes and Correspondence, 1988.
*/

#pragma once

#include "SaturationAdjustment.h"

#include <cmath>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <string>
#include <iostream>

class cNeptuneModel;
class Array;

using namespace std;

class SaturationAdjustmentNept {
public:

    explicit SaturationAdjustmentNept(cNeptuneModel& model) : m(model) {}

    // Selects between ATNEPT's inherited routine and the SHARED SaturationAdjustment<Planet>.
    // DEFAULT 1 = the shared algorithm SINCE 2026-10-10 (the user's decision, taken together with
    // the H2O and NH3 ice pairs in cNeptuneModel.h); ATNEPT_SATADJ=0 restores the inherited
    // routine. Every ATNEPT number before that commit was made with the inherited routine unless
    // its run set the knob, and with the liquid coefficients in the ice slots.
    //
    // 224 iterations, 8 threads (ATJUP/satchk/giants/NEPT224i, NEPT224s), inherited -> shared with
    // the ice pairs: 58 of 118 printed extrema identical, temperature, u and w among them;
    // max h2o_cloud 23.99 -> 24.73 g/m3, max nh3_ice 51.40 -> 61.17 (+19 %), max latent heat
    // 47.11 -> 45.75 W/m3; nh3_cloud and ch4_cloud, exactly zero with the inherited routine, reach
    // 1.43 and 5.45 g/m3. The shared routine's own budget over its 113 calls: no ice deleted,
    // no cell unconverged, column change of each gas below 1e-6 g/m2. The inherited routine has
    // no such instrument. Which of the two is right on Neptune is not settled by a measurement.
    static int mirrored_enabled(){
        static const int v = [](){
            const char* e = getenv("ATNEPT_SATADJ"); return e ? atoi(e) : 1; }();
        return v;
    }

    // Dispatch. The call-site signature is unchanged, so cNeptuneModel.cpp needs no edit.
    void run(const std::string& gas,
             double coeff_A,   double coeff_B,
             double coeff_A_i, double coeff_B_i,
             double t_0,       double t_00,
             double ep,        double lv,  double ls,
             double cp,        double r,
             double C,         double L0,  double R,
             double del_alf,   double del_bet,   double m_mol,
             Array& c,         Array& cloud,   Array& ice);

    // ATNEPT's own algorithm, unchanged and still the default path.
    void run_legacy(const std::string& gas,
             double coeff_A,   double coeff_B,
             double coeff_A_i, double coeff_B_i,
             double t_0,       double t_00,
             double ep,        double lv,  double ls,
             double cp,        double r,
             double C,         double L0,  double R,
             double del_alf,   double del_bet,   double m_mol,
             Array& c,         Array& cloud,   Array& ice);

    // -----------------------------------------------------------------------
    // Static helper functions — usable without a SaturationAdjustmentNept instance
    // -----------------------------------------------------------------------
    static double clausius_clapeyron(double T_K, double A, double B){
        return std::exp(A / T_K + B);
    }

    static double saturation_vapour_pressure(double T_K,
            double C, double L0, double R, double del_alf, double del_bet){
        return std::exp(C
            + (-L0 / T_K + del_alf * std::log(T_K) + del_bet * T_K)
            / (1e-3 * R));
    }

    static double humility_critical(double x, double Hu_cr_max, double Hu_cr_mid){
        return (Hu_cr_max - Hu_cr_mid) * (x * x - 2.0 * x) + Hu_cr_max;
    }

private:
    cNeptuneModel& m;

    static constexpr int iter_prec_end = 30;
    static constexpr double q_diff_min = 1.0e-4;
};
