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

    // Selects between ATNEPT's inherited routine and the SHARED SaturationAdjustment<Planet>
    // that ATSAT and ATJUP instantiate. Default 0 = the inherited one, so the model is
    // byte-identical until this is set. ATNEPT_SATADJ=1 selects the shared algorithm.
    //
    // NOT YET EVALUATED ON NEPTUNE. On ATSAT the switch was measured to condense 16 % less peak
    // cloud water and move the deck a layer higher, and which of the two is right was not settled
    // there either. The same caveat applies here and one more besides: the ice-phase coefficient
    // quadruple passed below is the LIQUID one, because ATNEPT's parameter set has no ice pair
    // for H2O, NH3 or CH4 — and on Neptune, colder than Saturn, the ice branch is the one that
    // matters most. Supplying real ice coefficients is what makes this knob worth turning on.
    static int mirrored_enabled(){
        static const int v = [](){
            const char* e = getenv("ATNEPT_SATADJ"); return e ? atoi(e) : 0; }();
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
