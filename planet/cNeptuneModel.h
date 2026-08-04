#ifndef _CNEPTUNEMODEL_H
#define _CNEPTUNEMODEL_H

#include <fenv.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Array.h"
#include "BoundaryConditions.h"   // BCForm, and the shared BC template
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"
#include "PythonStream.h"
#include "Utils.h"
#include "Config.h"

// Forward declarations for standalone friend classes
class BC_Nept;
class ChemistryNept;
class SaturationAdjustmentNept;
class PressureSolverNept;
class VelocityInitializerNept;


#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

namespace{
    std::function<double(double)> default_lambda=[](double i)->double{return i;};
}

class cNeptuneModel{

    // Shared physics/output templates (SHARED.md5, `make check-shared`). ATNEPT is the third
    // model to take these; see ParaViewWriter.h for what it provides and what it does not.
    template<class M> friend class ParaViewWriter;
    template<class M> friend class ConvectiveAdjustment;
    template<class M> friend class FluxLimiter;
    template<class M> friend class SaturationAdjustment;
    template<class M> friend class PressureSolver;
    template<class M> friend class Reporting;
    template<class M> friend class Radiation;
    template<class M> friend class BoundaryConditions;
    template<class M> friend class Turbulence;
    template<class M> friend class Precipitation;
    friend class BC_Nept;
    friend class ChemistryNept;
    friend class SaturationAdjustmentNept;
    friend class PressureSolverNept;
    friend class VelocityInitializerNept;

public:

    // ---- Hooks for the shared ParaViewWriter<Planet> (ParaViewWriter.h) ----
    // planet_name() is the word in an output FILE name ("Neptune_radial_20_1.vtk");
    // planet_short() is the abbreviation inside a .vtk title line
    // ("Radial_Data_Nept_Circulation"). ATNEPT carried both spellings by hand, and its
    // "has been written" line used the SHORT one — so it announced Nept_radial_20_1.vtk,
    // a file that does not exist. The shared writer names the file it actually wrote.
    // The model's own name in log lines written by SHARED code — "ATNEPT: ..." — so a
    // shared header can say which planet it is running on without knowing anything else
    // about it. ATSAT and ATJUP carry the same accessor.
    static const char* planet_tag(){ return "ATNEPT"; }

    /*
     * ---- Temperature bounds for the RK4 integrator, as PHYSICAL temperatures ----
     *
     * These were bare literals in RungeKutta_Nept_Turb.cpp — t_min = 0.1, t_max = 10.0 — in
     * NONDIMENSIONAL units, with comments reading "~7.6 K" and "~760 K physical". Those comments
     * were right for a t_ref they were written against and wrong here: Neptune's t_ref is 72.5, so
     * the same constants mean 7.25 K and 725 K.
     *
     * 725 K IS BELOW THIS MODEL'S OWN INITIAL PROFILE. init_temperature builds the deep equator at
     * 1111 K, so the first RK step clamped every interior cell there down to the ceiling while
     * i = 0 — outside the integrator's i = 1..im-2 range — kept its 1111 K. Measured: max|t-t_init|
     * went 0 -> 377.664 K in one step, at i = 1, and stayed at exactly that value on the next step,
     * which is what a clamp looks like and what a tendency does not.
     *
     * Expressed in KELVIN here and divided by t_ref at the point of use, so the number means the
     * same thing on every planet and cannot go stale when t_ref changes. The ceiling is a NUMERICAL
     * guard against a runaway, not a physical claim: it is set clear of anything the model
     * constructs rather than at a temperature Neptune is believed to reach.
     */
    static double t_min_K(){ return 7.5;    }   // floor; prevents a buoyancy blow-up
    static double t_max_K(){ return 2000.0; }   // ceiling; well above the 1111 K init_temperature builds

    // ---- What the SHARED BoundaryConditions.h asks of this model ----
    //
    // The field lists, the loop margin, the default extrapolation form and each knob's default.
    // Every one is a model FACT rather than a variant of the algorithm, which is why the shared
    // header asks rather than assumes. Two of ATNEPT's answers differ from ATSAT's and both are
    // deliberate:
    //
    //   bc_margin() = 0.  ATNEPT applies its boundary conditions over the FULL j,k ranges, where
    //   ATSAT works the interior rows only (margin 1). Changing that would change which cells the
    //   corners get, so it is stated rather than harmonised.
    //
    //   bc_default_form() = NEUMANN, the (4/3,-1/3) two-point extrapolation, where ATSAT defaults
    //   to the three-point CUBIC. This is not an oversight: BC_Nept.h records the reason — the
    //   cubic "amplifies alternating errors by 7x per call and is unstable near the SeaMount
    //   contour". ATJUP defaults to NEUMANN for its own reasons. BCForm::DEFAULT = 0 meaning
    //   "this planet's own form" is exactly what lets three models disagree here without the
    //   shared code choosing.
    //
    // Every hardening knob is OFF, because on Neptune none of them has been measured.
    static int bc_margin(){ return 0; }
    static int bc_default_form(){ return BCForm::NEUMANN; }
    static int bc_default_rigid_lid(){ return 0; }
    static int bc_default_top_taper(){ return 0; }
    static int bc_default_pole_copy(){ return 0; }
    static int bc_default_radius_copy(){ return 0; }

    // Snapshot of the lid temperature, for the shared BoundaryConditions' lid-pin knob. It is
    // DECLARED here and filled by the shared header itself, lazily, on the first call with the pin
    // enabled — which is why there is no initialisation to write: with the knob off (Neptune's
    // default, as on Saturn) it stays empty and costs one size() test per call.
    //
    // Worth knowing before anyone enables it: on ATSAT the pin was measured to be a cure for a
    // drift that does not exist — the lid moved +0.077 K over 28 iterations — and, because the
    // snapshot is taken during INITIALISATION, pinning holds a pre-first-iteration value rather
    // than "where the lid would otherwise have been". Neptune's lid has not been measured at all.
    // ---- Turbulence closure fields, for the SHARED Turbulence.h ----
    //
    // STAGE ONE OF THREE, following the staging ATSAT used (902c609, 4b1072c, then the coupling).
    // What arrives with this set is the CLOSURE: it reads the velocity field and fills nue* and
    // its diagnostics. k* and dis* are allocated and written by the closure but are NOT yet
    // prognostic — nothing integrates them in RungeKutta_Nept, which predates the closure
    // entirely — and nue* reaches no momentum or scalar equation. Those are stages two and three.
    //
    // tken/disn exist now so the prognostic stage has the start-of-step copies it will need
    // without a second pass over this header.
    // ---- RK4 stage accumulators, for the separated integrator ----
    //
    // The running sum k1 + 2k2 + 2k3 + k4 needs a place to live once the four stages stop sharing
    // one cell loop. RungeKutta_Nept currently runs all four stages inside a single parallel loop
    // over cells, holding k1..k4 as per-thread scalars, while RHSNept differentiates t,u,v,w at
    // i+-1, j+-1, k+-1 — cells other threads are simultaneously overwriting. That is why two runs
    // of the same binary at 24 threads differ in 4 of 7 output files, measured.
    //
    // These are declared and allocated now, inert, so that the integrator rewrite is a change to
    // one file. ATSAT fixed the same defect the same way in 71082e7 and became reproducible at any
    // thread count; that is the acceptance test this is aiming at.
    Array acc_tke;              // RK4 accumulator for k*
    Array acc_dis;              // RK4 accumulator for dis*
    Array acc_t;
    Array acc_u;
    Array acc_v;
    Array acc_w;
    Array acc_ch4;
    Array acc_ch4_cloud;
    Array acc_ch4_ice;
    Array acc_h2o;
    Array acc_h2o_cloud;
    Array acc_h2o_ice;
    Array acc_h2s;
    Array acc_h2s_cloud;
    Array acc_h2s_ice;
    Array acc_nh3;
    Array acc_nh3_cloud;
    Array acc_nh3_ice;
    Array acc_nh4sh;

    Array P_rain;
    Array P_snow;
    Array P_graupel;
    Array P_nh3_rain;
    Array P_nh3_snow;
    Array P_nh3_graupel;
    Array P_ch4_rain;
    Array P_ch4_snow;
    Array P_ch4_graupel;
    Array P_nh4sh;
    Array Q_precip;
    Array S_precip_h2o;
    Array S_precip_h2o_cloud;
    Array S_precip_h2o_ice;
    Array S_precip_nh3;
    Array S_precip_nh3_cloud;
    Array S_precip_nh3_ice;
    Array S_precip_ch4;
    Array S_precip_ch4_cloud;
    Array S_precip_ch4_ice;
    Array_2D precip_srf_h2o;
    Array_2D precip_srf_nh3;
    Array_2D precip_srf_ch4;
    Array_2D precip_srf_nh4sh;
    Array_2D precip_srf_total;
    Array rhs_tke;              // tendency of k*,   assembled in RHS_Nept
    Array rhs_dis;              // tendency of dis*, assembled in RHS_Nept
    // VERIFIED TO EVOLVE, and three commits claimed otherwise. With ATNEPT_TURB=1 the maximum of
    // k* grows monotonically over six iterations — 0.000779, 0.000811, 0.000841, 0.000871,
    // 0.000900, 0.000928 — and a direct probe shows the expected mixture of signs in rhs_tke
    // (-1.29e-03 at i=1, +1.24e-03 at i=2). Stage two has worked since it landed.
    //
    // f7c3e96, 07252d9 and aac08d7 said it did not, on a measurement error of mine: the report
    // prints max AND min on one line, and `sed 's/.*= *//'` is greedy, so it returned the MIN
    // column — which is zero at the poles. Any extraction from this report must say which column
    // it wants.
    //
    // One real observation survives from that detour: k* holds its k_abl value in the two layers
    // below abl_height = 20 km (7.470e-04) and sits on the background floor above (7.649e-05), so
    // only 2 of Neptune's 41 levels are inside the boundary layer. abl_height is 20000 m in all
    // three models, which is worth revisiting for grids this different — but it is not a defect
    // and nothing failed because of it.
    Array tke;                  // turbulent kinetic energy k*      [dimensionless]
    Array dis;                  // dissipation eps* or omega*       [dimensionless]
    Array tken;                 // k* at the start of the RK4 step
    Array disn;                 // dis* at the start of the RK4 step
    Array nue;                  // eddy viscosity nue* (the closure's own name)
    Array nue_t;                // eddy viscosity as an RHS would read it
    Array prod;                 // shear production P_k
    Array tke_source;           // P_k - Y_k
    Array dis_source;           // P_w - Y_w + D_w
    Array_2D vel_star;          // friction velocity u_tau at the first fluid layer [m/s]

    double re_turb = 1.0;       // = vel_star_ref*z_0/nue, set by the closure
    double abl_height = 20000.0; // boundary-layer height [m]

    // turb_model is a configuration parameter, as in ATJUP, ATSAT and ATURAN — declared by
    // NeptuneParams.h.inc from param.py, so it is NOT declared here. ATNEPT_TURB_MODEL still
    // overrides it at runtime.
    bool turb_active = false;   // THE gate; set in Run() from ATNEPT_TURB and turb_model together

    std::vector<std::vector<double> > t_top_init;

    std::vector<Array*> bc_fields_radius();
    std::vector<Array*> bc_fields_theta_extrap();
    std::vector<Array*> bc_fields_theta_zero();
    std::vector<Array*> bc_fields_phi();

    // ATNEPT has no turbulence fields at all — Turbulence.h is not among its shared headers and
    // its RungeKutta predates the closure — so the turbulence boundary pass has nothing to act on
    // and is switched off at the source rather than given empty work.
    std::vector<Array*> bc_turb_fields(){ return {}; }
    std::vector<double> bc_turb_floors(){ return {}; }
    bool bc_turb_active() const { return false; }

    /*
     * ---- Neptune's radiative constants, for the SHARED Radiation.h ----
     *
     * MEASURED PROPERTIES OF NEPTUNE, which is why they live here and not in the shared file: a
     * mechanical copy of Saturn's radiation would have produced Saturn's budget on Neptune's grid.
     *
     *                          ATJUP    ATSAT    ATNEPT   source
     *   solar constant         50.5     14.83    1.505    1361/a^2, a = 5.20, 9.58, 30.07 AU
     *   Bond albedo            0.343    0.342    0.290    Pearl & Conrath
     *   intrinsic flux F_int   5.4      2.01     0.433    Pearl & Conrath 1991, 0.433 +/- 0.046
     *   H2 mole fraction       0.86     0.96     0.80
     *   He mole fraction       0.136    0.032    0.19
     *
     * Neptune is the extreme case of the set: it receives ~1/10 of Saturn's sunlight and ~1/34 of
     * Jupiter's, yet radiates about 2.6x what it absorbs, so its INTERNAL flux dominates its budget
     * far more than either. Whether this scheme, calibrated on Jupiter, behaves sensibly in that
     * regime is exactly what the knob is for and is not established here.
     */
    static double rad_F_int()       { return 0.433; }  // Neptune intrinsic heat flux [W/m2]
    static double rad_S_solar()     { return 1.505; }  // solar constant at 30.07 AU [W/m2]
    static double rad_albedo_bond() { return 0.290; }  // Neptune Bond albedo (S*(1-A) is ABSORBED)
    static double rad_x_H2()        { return 0.80;  }  // H2 mole fraction
    static double rad_x_He()        { return 0.19;  }  // He mole fraction

    // Physical thickness of layer i in metres, which the radiation integrates optical depth over.
    double layer_thickness_m(int i){
        if(i < 0 || i > im-2) return 0.0;
        return (double)(m_layer_heights[i+1] - m_layer_heights[i]) * 1.0e3;
    }

    // Neptune has no surface: the column starts at i = 0 everywhere.
    int surface_index(int, int) const { return 0; }

    // p_dyn is stored as the NONDIMENSIONAL kinematic pressure, so displaying it in bar needs
    // r_mix*u_0^2*1e-5. Default OFF (returns 1.0) so the printed number is unchanged; ATNEPT_PDYN_UNITS=1
    // makes it actually bar. Mirrors ATJUP's accessor of the same name, which is what lets the
    // pressure row read the same in all three models.
    double p_dyn_to_bar() const {
        static const bool on = [](){ const char* e = getenv("ATNEPT_PDYN_UNITS"); return e && atoi(e) != 0; }();
        return on ? r_mix * u_0 * u_0 * 1.0e-5 : 1.0;
    }


    // ---- ATSAT_QHEAT_SCALE / the latent+sensible heating scale ----
    //
    // Returns exactly 1.0 when off, so the default path is bit-identical.
    //
    // Thermo_*.cpp forms the heating as  lv * velocity_av * grad(rho) / (L_atm * L_atm), where
    // velocity_av is built from the RAW NONDIMENSIONAL u,v,w and grad() is per nondimensional
    // length. The physical volumetric rate is lv * (v . grad rho) in W/m3, which needs u_0 to
    // make the velocity a speed and ONE division by the length scale IN METRES. The code divides
    // by the length scale twice, in kilometres. The ratio is u_0 * L_atm / 1e3.
    //
    // Q_Latent and Q_Sensible are OUTPUT-ONLY — nothing in RHS_*/RungeKutta_* reads them — so
    // this changes no physics, but it does change the .vtk/.vts values, which is why it is a knob
    // and not a silent repair.
    double qheat_fix() const {
        static const int on = [](){ const char* e = getenv("ATNEPT_QHEAT_SCALE"); return e ? atoi(e) : 0; }();
        return on ? u_0 * L_atm / 1.0e3 : 1.0;
    }


    // ---- Hooks for the shared FluxLimiter<Planet> (FluxLimiter.h) ----
    //
    // metricRadius() is the established hook for the one place ATSAT and ATJUP genuinely differ
    // in the limiter: ATJUP shifts rad.z itself at initialisation and so returns rm unchanged,
    // while ATSAT shifts the metric factors here instead. ATNEPT is in ATJUP's position for a
    // simpler reason — it has no metric radius at all — so this is the identity unless
    // ATNEPT_METRIC_RADIUS is set, and the shared limiter reproduces the m.rad.z[i] the
    // hand-written copy used, exactly.
    double metricRadius(double rm){
        static const double R_km = [](){
            const char* e = getenv("ATNEPT_METRIC_RADIUS"); return e ? atof(e) : 0.0; }();
        if(!(R_km > 0.0)) return rm;
        return rm + (R_km / L_atm - 1.0);
    }

    // ---- Hooks for the shared Reporting<Planet> (Reporting.h) ----

    // ATJUP flips cos(theta) in the southern hemisphere for the continuity residual; ATSAT and
    // ATNEPT never have. Same knob shape as ATJUP's, so the three models answer one question
    // rather than differ by a missing line. Default off = unchanged.
    static bool costhe_abs(){
        static const bool v = [](){ const char* e = getenv("ATNEPT_COSTHE_ABS"); return e && atoi(e) != 0; }();
        return v;
    }

    // Column layout of the min/max report. ATNEPT uses the same widths ATSAT does — 6 for the
    // unit column, ten spaces between the max and min halves. ATJUP widened its unit column to 12
    // and uses three spaces, because its unit strings are longer.
    static int minmax_unit_width()      { return 6; }
    static const char *minmax_separator(){ return "          "; }

    static const char *steady_heading(){
        return " 3D iterational process for the surface boundary conditions\n printout of maximum and minimum absolute and relative errors of the computed values at their locations: level, latitude, longitude";
    }

    // The iteration line of the steady-state header, including its trailing newline.
    //
    // THIS PRINTS iter_n, NOT n, AND THAT IS A CORRECTION. ATNEPT's steadyQuery printed `n`, a
    // member that is DECLARED AND NEVER ASSIGNED ANYWHERE in this model — the iteration loop
    // counts with iter_n. It therefore printed stack garbage (observed: "n = 2185456" on a run
    // whose nm was 2). Nobody had seen it because nothing called the routine. ATJUP genuinely
    // counts with n and prints it; ATSAT and now ATNEPT use iter_n.
    std::string steady_iter_line() const {
        return "      n = " + std::to_string(iter_n) + "\n";
    }

    // ---- What the SHARED PressureSolver.h asks of this model ----
    //
    // has_obstacle() is the same fact is_solid() states cell by cell, asked once: Neptune contains
    // no solid body, so the solver can skip its obstacle handling entirely.
    static bool has_obstacle(){ return false; }

    // Whether the projection treats the radial walls as a rigid lid — u = 0 there rather than
    // extrapolated. False keeps ATNEPT's existing open boundaries. ATSAT answers false too; the
    // knob that turns it on there is a separate, still-unsettled question.
    static bool press_rigid_lid(){ return false; }

    // ATNEPT does not stretch the radial coordinate, so the solver's exp_rm factor stays 1.
    // ATJUP's coord_stretching forms 1/(rm+1), which only makes sense while rad.z starts at 1.
    bool coord_stretching = false;

    // Radial-wall conditions on the intermediate velocity and the RHS, applied before the
    // pressure Poisson solve. Mirrored from ATSAT's, which is where the cubic (4/3,-1/3)-style
    // extrapolation and the rigid-lid alternative are explained. Defined in Pressure_Nept.cpp.
    void prepareProjectionBoundaries(bool rigid_lid);

    // ---- What the SHARED SaturationAdjustment.h asks of this model ----

    // Neptune contains no solid obstacle: it is a gas giant modelled as a spherical shell, and
    // nothing in ATNEPT marks a cell as ground. ATSAT answers the same question the same way;
    // only ATJUP has an obstacle, and only for its seamount experiments.
    bool is_solid(int, int, int) const { return false; }

    // The density the adjustment divides by. ATNEPT_LOCAL_RHO=1 uses the LOCAL mixture density
    // where it is usable and falls back to the constant r_mix where it is not — which matters
    // because rho_mix is zero until computeMixtureDensity has run, and a zero here would divide
    // through the whole adjustment. Default 0 = the constant r_mix everywhere, which is what
    // ATNEPT has always used. Exactly ATSAT's and ATJUP's rho_at().
    double rho_at(int i, int j, int k){
        static const int local = [](){
            const char* e = getenv("ATNEPT_LOCAL_RHO"); return e ? atoi(e) : 0; }();
        if(!local) return r_mix;
        const double rho = rho_mix.x[i][j][k];
        return (rho > 0.0 && std::isfinite(rho)) ? rho : r_mix;
    }

    // Whether the shared adjustment writes the static pressure back after condensing. ATSAT says
    // true; ATNEPT's own routine does not touch p_stat, so it says false and the shared algorithm
    // leaves the field alone — one of the two models' behaviours had to be named rather than
    // assumed.
    static bool satadj_updates_pstat(){ return false; }

    // The model's own floor on sin(theta) in the METRIC. ATNEPT declares none, as ATSAT does not;
    // it exists so ATPhys::polar_divisor_floor<Planet>() compiles. With ATNEPT_SINTHE_TRACK unset
    // that function returns the literal 0.4 the hand-written limiter used, so this value is not
    // reached by default. See the ATSAT commit that made that floor one accessor instead of four
    // drifting literals.
    static double sinthe_min(){
        static const double v = [](){
            const char* e = getenv("ATNEPT_SINTHE_MIN");
            const double x = e ? atof(e) : 0.0;
            return (x >= 0.0 && x < 1.0) ? x : 0.0;
        }();
        return v;
    }

    static const char* planet_name(){ return "Neptune"; }
    static const char* planet_short(){ return "Nept"; }

    // What the panorama .vts prints in its "Temperature" array. Neptune writes degrees
    // Celsius, as ATJUP does; ATSAT writes kelvin/10. See the note on the same hook in
    // cSaturnModel.h — one array name, three models, two different quantities.
    double paraview_temperature(double t_nd) const { return t_nd * t_ref - 273.15; }


    const char *filename;

    cNeptuneModel();
    ~cNeptuneModel();

    cNeptuneModel(const cNeptuneModel&) = delete;
    cNeptuneModel& operator=(const cNeptuneModel&) = delete;

    static cNeptuneModel* get_model(){
        if(!m_model){
            m_model = new cNeptuneModel();
        }
        return m_model;
    }

    void LoadConfig(const char *filename);
    void Run();

    #include "NeptuneParams.h.inc"

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0;


    int n, n_print, n_paraview, panorama, panorama_step;
    int iter, iter_max, j_max;
    int velocity_iter, pressure_iter;
    int iter_n, panorama_cnt;
    int i_res, j_res, k_res;

    std::vector<double> tropopause_layers; // keep the tropopause layer index
    std::vector<std::vector<int> > i_topography;
    std::vector<double> u_trans;
    std::vector<double> v_trans;
    std::vector<double> w_trans;

    double maxValue, minValue;
    /*
     * Given a latitude, return the layer index of tropopause
    */
    int get_tropopause_layer(int j){
        assert(j>=0);
        assert(j<jm);
        //refer to  BC_Thermo::TropopauseLocation and BC_Thermo::GetTropopauseHightAdd
        //tropopause height is proportional to the mean tropospheric temperature.
        //higher near the equator - warm troposphere
        //lower at the poles - cold troposphere

//        tropopause_layers[j] = 30;
//        tropopause_layers[j] = 35;
//        init_tropopause_layers();
        tropopause_layers[j] = im_tropopause[j];
        return tropopause_layers[j];
    }
    /*
     *
    */
    int get_surface_layer(int j, int k){
        return i_topography[j][k];
    }
    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    std::vector<float> get_layer_heights(){
        return m_layer_heights;
    }   
    /*
    * Given a altitude, return the layer index
    */
    int get_layer_index(float height){
        std::size_t i = 0;
        for(; i<m_layer_heights.size(); i++){
            if(height<m_layer_heights[i])
                return i-1;
        }
        return i;
    }



private:

    static cNeptuneModel* m_model;

    PythonStream ps;
    std::streambuf *backup;

    const double c43 = 4.0/3.0, c13 = 1.0/3.0, c32 = 3.0/2.0, c42 = 4.0/2.0, c12 = 1.0/2.0;
    static const int im = 41, jm = 181, km = 361;

    int i_max = 40;  // corresponds to about 125 km above 10e6 Pa pressure level, maximum hight of the tropopause at equator
    int i_beg = 30;  // corresponds to about 100 km above 10e6 Pa pressure level, maximum hight of the tropopause at poles
//    int i_beg = 40;  // corresponds to about 100 km above 10e6 Pa pressure level, maximum hight of the tropopause at poles

    double gam = 0.0;    // dry adiabatic lapse rate [K/km], computed as g*1e3/cp_mix after ThermalPropertiesUran()
    double mue_mix, k_mix, cp_mix, r_mix, R_mix, c_mix, x_mix;

// at 230K NH3 and H2S condense via a heterogenious reaction: NH3 + H2S -> NH4SH ( Planetary Scienses, de Pater, Lissauer)


// temperatures at triple point and ice formation
    double t_0_h2o = 273.15;  // in K == 0°C, triple point
    double t_00_h2o = 210.15;  // in K == -67°C (Planetary Siences)
//    double t_00_h2o = 241.15;  // in K == -32°C (COSMO)
    double t_000 = 235.15;  // in K == -20°C (precipitation module)

    double t_0_h2s = 187.65;  // in K == -85.5°C, h2s-ice cloud formation (Planetary Sciences, p. 96)
    double t_00_h2s = 220.0;  // in K == -133.15°C, h2s-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_ch4 = 90.69;  // in K == -182.456°C, triple point
    double t_00_ch4 = 190.56;  // in K == -85.5°C, h2s-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh3 = 195.5;  // in K == -77.65°C, triple point, gas and liquid pase
    double t_00_nh3 = 220.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)
//    double t_00_nh3 = 140.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)
//    double t_00_nh3 = 120.0;  // in K == -153.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)

//    double t_0_nh4sh = 230.0;  // in K == -43.15°C, nh4sh formation onset (Planetary Sciences, p. 96)
//    double t_00_nh4sh = 200.0;  // in K == -73.15°C, nh4sh formation end (Planetary Sciences, p. 96)

    double t_0_nh4sh = 230.0;  // in K == -57.15°C, nh4sh formation onset (Planetary Sciences, p. 96)
    double t_00_nh4sh = 200.0;  // in K == -113.15°C, nh4sh formation end (Planetary Sciences, p. 96)


// pressures take from  Planetary Sciences, p. 96
    double p_0_h2o = 21.0;  // in bar
    double p_00_h2o = 50.0;  // in bar

    double p_0_ch4 = 0.1;  // in bar
    double p_00_ch4 = 1.1;  // in bar

    double p_0_h2s = 2.8;  // in bar
    double p_00_h2s = 7.0;  // in bar

    double p_0_nh3 = 20.2;  // in bar
    double p_00_nh3 = 38.0;  // in bar

    double p_0_nh4sh = 20.5;  // in bar
    double p_00_nh4sh = 38.0;  // in bar


// constants for Clausius-Clapeyron law
    double coeff_h2_A = -2000.0; // invented
    double coeff_h2_B = 8.0; // invented

    double coeff_he_A = -3000.0; // invented
    double coeff_he_B = 10.0; // invented


    double coeff_h2o_A = -4961.04;  // from triple and critical point values for water
    double coeff_h2o_B = 13.0662;  // from triple and critical point values for water

    double coeff_h2o_A_i = -4961.04; // invented
    double coeff_h2o_B_i = 13.0662; // invented


    double coeff_ch4_A = -1033.3;  // from triple and critical point values for methane
    double coeff_ch4_B = 6.3910;  // from triple and critical point values for methane

    double coeff_ch4_A_i = -1033.3; // invented
    double coeff_ch4_B_i = 6.3910; // invented

    double coeff_h2s_A = -2251.66;  // from triple and critical point values for hydrogen sulfide
    double coeff_h2s_B = 10.5253;  // from triple and critical point values for hydrogen sulfide

    double coeff_h2s_A_i = -2251.66; // invented
    double coeff_h2s_B_i = 10.5253; // invented


    double coeff_nh3_A = -2836.56;  // from triple and critical point values for ammonia
    double coeff_nh3_B = 11.7271;  // from triple and critical point values for ammonia

    double coeff_nh3_A_i = -2836.56; // invented
    double coeff_nh3_B_i = 11.7271; // invented


    double coeff_nh4sh_A = -2836.56;  // invented
    double coeff_nh4sh_B = 11.7271; // invented


// condensate/crystal densities in kg/m³ (Planetary Sciences p. 90, 2010)
// These are phase-change densities for cloud microphysics — NOT for gas-phase diffusion.
    double rho_cond_h2    = 0.408;    // liquid hydrogen density in kg/m³
    double rho_cond_he    = 0.1752;   // liquid helium density in kg/m³
    double rho_cond_nh3   = 0.7623;   // liquid ammonia density in kg/m³
    double rho_cond_nh4sh = 1170.0;   // solid ammonium hydrosulfide density in kg/m³
    double rho_cond_h2s   = 1.5357;   // liquid hydrogen sulfide density in kg/m³
    double rho_cond_h2o   = 1000.0;   // liquid water density in kg/m³
    double rho_cond_ch4   = 0.657;    // liquid methane density in kg/m³
    double rho_cond_mix   = 0.0;      // computed in ThermalPropertiesNept

// mass density (concentration) of species in g/cm³ == 10e6 g/m³ == 10e3 kg/m³                            Planetary sciences p. 90 2010
    double rg_h2 = 0.408; // mass density of vapour
    double rg_he = 0.1752; // mass density of vapour
    double rg_nh3 = 0.7623; // mass density of vapour
    double rg_nh4sh = 1170.0; // mass density of vapour
    double rg_h2s = 1.5357; // mass density of vapour
    double rg_h2o = 0.005; // mass density of vapour
    double rg_ch4 = 0.657; // mass density of vapour
    double rg_mix = 0.0;   // computed in ThermalPropertiesNept

// NH4SH crystal radius for Stokes settling
    double r_p_nh4sh = 1.0e-5;  // NH4SH crystal radius in m

// molecular weights
    double m_h2 = 2.016;  // molecular weight of hydrogen in kg/Kmol (molar mass)
    double m_he = 4.02602;  // molecular weight of helium in kg/Kmol
    double m_nh3 = 17.03052;  // molecular weight of ammonia in kg/Kmol
    // Reaction enthalpy of NH3(g) + H2S(g) -> NH4SH(s), per kg of NH4SH formed.
    //
    // From standard enthalpies of formation:
    //     dHf NH3(g)    -45.90 kJ/mol
    //     dHf H2S(g)    -20.60 kJ/mol
    //     dHf NH4SH(s) -156.90 kJ/mol
    //     dH_rxn = -156.90 - (-45.90 - 20.60) = -90.40 kJ/mol   (exothermic)
    // divided by m_nh4sh = 51.1114 kg/kmol:  90.40e3 / 0.0511114 = 1.7687e6 J/kg.
    //
    // Sign convention here is HEAT RELEASED, so the source term is + dh_nh4sh * w_nh4sh and is
    // positive where NH4SH is forming. w_nh4sh is kg/(m3 s), so the product is W/m3 — the same
    // dimension as the enthalpy-diffusion term it is added to, which is the whole point.
    double dh_nh4sh = 1.7687e6;   // reaction enthalpy of NH4SH formation in J/kg
    double m_nh4sh = 51.1114;  // molecular weight of ammonium hydrosulfide in kg/Kmol
    double m_h2s = 34.08088;  // molecular weight of hydrogen sulfide in kg/Kmol
    double m_h2o = 18.01588;  // molecular weight of water in kg/Kmol
    double m_ch4 = 16.042;  // molecular weight of water in kg/Kmol

// vapour mass densities of gases                Planetary sciences p. 90 2010
    double r_h2 = 0.864;  // density of hydrogen vapour in kg/m³
    double r_he = 0.136;  // density of helium vapour  in kg/m³
    double r_ch4 = 0.19;  // density of ammonia vapour in kg/m³
    double r_nh3 = 0.1;  // density of ammonia vapour in kg/m³
//    double r_h2s = 0.8;  // density of hydrogen sulfide vapour in kg/m³
    double r_h2s = 0.3;  // density of hydrogen sulfide vapour in kg/m³
    double r_h2o = 0.08;  // density of water vapour in kg/m³
    double r_nh4sh = 0.007;  // density of ammonium hydrosulfide vapour in kg/m³  assumption
    double r_nh3_add = 0.09;  // density of ammonia vapour in kg/m³

// vapour mass densities of clouds and ices               Planetary sciences p. 90 2010
    double r_nh3_cloud = 0.004;  // density of ammonia cloud in kg/m³
    double r_h2o_cloud = 0.009;  // density of water cloud in kg/m³
    double r_nh3_ice = 0.0003;  // density of ammonia ice in kg/m³
    double r_h2o_ice = 0.0038;  // density of water ice in kg/m³
    double r_h2s_ice = 0.007;  // density of water ice in kg/m³
    double r_ch4_ice = 0.0082;  // density of water ice in kg/m³

// vapour molar densities of gases
    double c_h2 = r_h2/m_h2;  // density of hydrogen vapour in kmol/m³
    double c_he = r_he/m_he;  // density of helium vapour  in kmol/m³
    double c_nh3 = r_nh3/m_nh3;  // density of ammonia vapour in kmol/m³
    double c_h2s = r_h2s/m_h2s;  // density of hydrogen sulfide vapour in kmol/m³
    double c_h2o = r_h2o/m_h2o;  // density of water vapour in kmol/m³
    double c_ch4 = r_ch4/m_ch4;  // density of water vapour in kmol/m³
    double c_nh4sh = r_nh4sh/m_nh4sh;  // density of ammonium hydrosulfide vapour in kmol/m³  assumption

 // ratio of vapour molecular weight to mean molecular weight
    double X_h2 = 0.864;
    double X_he = 0.136;
//    double X_h2o = 5.0e-5;
    double X_h2o = 1.7e-3;
    double X_nh3 = 2.0e-4;
    double X_h2s = 7.7e-5;
    double X_nh4sh = 3.6e-5;
    double X_ch4 = 3.6e-5;

// gas constants
    double R_h2 = 4124.2; // gas constant of hydrogen in J/(kg*K)
    double R_he = 2077.1; // gas constant of helium in J/(kg*K)
    double R_nh3 = 488.21; // gas constant of ammoinia in J/(kg*K)
    double R_nh4sh = 261.0; // gas constant of ammoinium hydrosulfide in J/(kg*K)   invented for the initial distribution of nh4sh
    double R_h2s = 243.96; // gas constant of hydrogen sulfide inJ/(kg*K) 
    double R_h2o = 461.52; // gas constant of water inJ/(kg*K)
    double R_ch4 = 518.28; // gas constant of methane inJ/(kg*K)

// dynamic viscosities
    double mue_h2 = 0.84e-5; // dynamic viscosity of hydrogen in Ns/m²
    double mue_he = 1.87e-5; // dynamic viscosity of helium in Ns/m²
    double mue_nh3 = 0.92e-5; // dynamic viscosity of ammonia in Ns/m²
    double mue_nh4sh = 0.99e-5; // dynamic viscosity of ammonium sulfide in Ns/m²
    double mue_h2s = 1.3e-5; // dynamic viscosity of hydrogen sulfid in Ns/m²
    double mue_h2o = 1.308e-3; // dynamic viscosity of water in Ns/m²
    double mue_ch4 = 1.107e-2; // dynamic viscosity of methane in Ns/m²

// thermal conductivities
    double k_h2 = 0.1317; // thermal conductivity of hydrogen in W/(m*K)
    double k_he = 0.1193; // thermal conductivity of helium in W/(m*K)
    double k_nh3 = 0.02102; // thermal conductivity of ammonia in W/(m*K)
    double k_nh4sh = 0.02102; // thermal conductivity of ammonium hydrosulfide in W/(m*K)
    double k_h2s = 0.013; // thermal conductivity of hydrogen sulfide in W/(m*K)
    double k_h2o = 0.0187; // thermal conductivity of water in W/(m*K)
    double k_ch4 = 0.0; // thermal conductivity of methane in W/(m*K)

// specific heat capacities
    double cp_h2 = 14.32e3;  // specific heat capacity of hydrogen in J/(kg*K)
    double cp_he = 5.19e3;  // specific heat capacity of helium in J/(kg*K)
    double cp_nh3 = 2.19e3;  // specific heat capacity of ammonia in J/(kg*K)
    double cp_nh4sh = 2.00e3;  // specific heat capacity of ammonium hydrosulfide in J/(kg*K)
    double cp_h2s = 2.24e3;  // specific heat capacity of hydrogen sulfid in J/(kg*K)
    double cp_h2o = 1.93e3;  // specific heat capacity of water in J/(kg*K)
    double cp_ch4 = 2.232e3;  // specific heat capacity of methane in J/(kg*K)

// ratios of gas constants of dry gas to vapour or vapour molecular weight to mean atmospheric molecular weight or m/m_mix
    double ep_h2 = 0.8572;  // ratio of the gas constants of dry air to h2 non-dimensional or m/m_mix
    double ep_he = 1.7152;  // ratio of the gas constants of dry he to h2 non-dimensional
    double ep_h2o = 8.1253;  // ratio of the gas constants of dry air to h2 non-dimensional
    double ep_h2s = 14.5192;  // ratio of the gas constants of dry hydrogen sulfide to h2 non-dimensional
    double ep_nh3 = 7.6752;  // ratio of the gas constants of dry ammonia to h2 non-dimensional
    double ep_ch4 = 7.6752;  // ratio of the gas constants of dry methane to h2 non-dimensional
    double ep_nh4sh = 21.7745;  // ratio of the gas constants of dry methane to h2 non-dimensional       invented

// latent heat of evaporation
    double lv_h2o = 2.5009e6;  // latent heat of h2o evaporation at 0°C in J/kg
    double lv_h2s = 3.5340e6;  // latent heat of h2s evaporation at -73°C in J/kg
    double lv_nh3 = 1.3720e6;  // latent heat of nh3 evaporation at -33.33 in J/kg
    double lv_ch4 = 5.11e5;  // latent heat of nh3 evaporation at -33.33 in J/kg

// latent heat of sublimation
    double ls_h2o = 2.8339e6;  // latent heat of h2o sublimation at 0°C in J/kg
    double ls_h2s = 7.4500e5;  // latent heat of h2s sublimation at -98°C in J/kg
    double ls_nh3 = 1.8320e6;  // latent heat of nh3 sublimation at -93.15 in J/kg
    double ls_ch4 = 5.11e5;  // latent heat of ch4 sublimation at -93.15 in J/kg

// Schmidt number
    double sc_h2 = 0.20;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_he = 0.22;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_ch4 = 0.99;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2o = 0.61;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2s = 0.94;  // Schmidt number of h2s, Sc = nue/D 
    double sc_nh3 = 0.61;  // Schmidt number of nh3, Sc = nue/D 
    double sc_nh4sh = 0.7;  // Schmidt number of nh4sh, Sc = nue/D 

// Sponge layer (top quarter of domain, quadratic Rayleigh damping on u)
    double alpha_sponge = 10.0;  // peak Rayleigh damping rate at top boundary

// Prantl numbers
    double Pr = 0.72;  // Prandtl-number 

// diffusion coefficients
    double D_nh3 = 1.5e-9; // ordinary diffusion coefficient of ammonia in m*m/s 
    double D_h2s = 1.36e-9; // ordinary diffusion coefficient of hydrogen sulfid in m*m/s 
    double D_nh4sh = 1.45e-9; // ordinary diffusion coefficient of ammonium hydrosulfide in m*m/s

// thermal diffusion coefficients
    double DT_nh3 = 1.54e-9; // thermal diffusion coefficient of ammonia in kg/(s*m)                           unklar
    double DT_h2s = 1.36e-9; // thermal diffusion coefficient of hydrogen sulfid in kg/(s*m)
    double DT_nh4sh = 1.45e-9; // thermal diffusion coefficient of ammonium hydrosulfide in kg/(s*m)

// constants for saturation vapour pressure and latent heat
    double C_h2o = 25.096;  //  in bar
    double C_ch4 = 1.627;  //  in bar
    double C_nh3 = 27.863;  //  in bar
    double C_h2s = 17.064;  //  in bar
    double C_nh4sh = 75.678;  //  in bar

    double L0_h2o = 3148.2;  //  in J/g
    double L0_ch4 = 553.1;  //  in J/g
    double L0_nh3 = 2016.0;  //  in J/g
    double L0_h2s = 747.0;  //  in J/g
    double L0_nh4sh = 2915.7;  //  in J/g

    double del_alf_h2o = 0.0;
    double del_bet_h2o = - 8.7e-3;

    double del_alf_ch4 = 1.002;
    double del_bet_ch4 = -4.1e-3;

    double del_alf_nh3 = - 0.888;
    double del_bet_nh3 = 0.0;

    // Ice-phase saturation quadruples, for the SHARED Precipitation.h. THESE MUST SIT AFTER THE
    // LIQUID CONSTANTS THEY COPY: member initialisers run in DECLARATION order, so declaring them
    // earlier reads uninitialised memory.
    //
    // CH4's are real values; H2O's and NH3's are the LIQUID constants standing in, exactly as
    // ATSAT's are — ATNEPT's parameter set has no ice pair for either, and inventing numbers for
    // Neptune's ices is a physics decision, not a port.
    double C_ch4_ice = 1.627;
    double del_alf_ch4_ice = 1.002;
    double del_bet_ch4_ice = -4.1e-3;
    double L0_ch4_ice = 553.1;
    double C_h2o_ice       = C_h2o;
    double del_alf_h2o_ice = del_alf_h2o;
    double del_bet_h2o_ice = del_bet_h2o;
    double L0_h2o_ice      = L0_h2o;
    double C_nh3_ice       = C_nh3;
    double del_alf_nh3_ice = del_alf_nh3;
    double del_bet_nh3_ice = del_bet_nh3;
    double L0_nh3_ice      = L0_nh3;

    double del_alf_h2s = 0.0;
    double del_bet_h2s = - 2.9e-3;

    double del_alf_nh4sh = - 1.760;
    double del_bet_nh4sh = 7.8e-4;
 
    std::vector<std::vector<int> > j_ellipse;
    bool has_welcome_msg_printed;
    double out_maxValue() const;
    double out_minValue() const;

    void init_layer_heights(){
        float h = L_atm/(im-1);
        for(int i=0; i<im; i++){
            m_layer_heights.push_back(i * h);
        } 
        return;
    }



    struct CellGeometry {
        // exp_rm/exp_2_rm are the radial coordinate-stretching factors the SHARED
        // PressureSolver.h reads. ATNEPT does not stretch (coord_stretching = false), so the
        // solver sets both to 1 and they are carried only so the struct satisfies the template.
        double rm, rm2, exp_rm, exp_2_rm;
        double sinthe, sinthe2, costhe, cotanthe;
        double inv_rm, inv_rm2;
        double inv_rmsinthe, inv_rm2sinthe, inv_rm2sinthe2;
        double costhe_inv_rm2sinthe;
        double inv_2dr, inv_2dthe, inv_2dphi;
        double inv_dr2, inv_dthe2, inv_dphi2;
    };

    void SetDefaultConfig();
    void RHSNept(int i, int j, int k, const CellGeometry& geo);
    void RungeKuttaNept();
    void NeptunePlotData();
    void paraview_vtk_longal(int n, int j_longal);
    void paraview_vtk_radial(int n, int i_radial);
    void paraview_vtk_zonal(int n, int k_zonal);
    void paraview_panorama_vts(int n);
    void paraview_sphere_vts(int n);

    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1., 
        std::function< double(double) > lambda = default_lambda,
        bool print_heading=false);

    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.0);

    void print_welcome_msg();
    void print_final_msg();
    void printMinMax();
    void initMsg();
    void writeResults();
    void writeData();

    void resetArrays();
    void TropopauseLocation();
    void Neptune_PlotData();


    void init_temperature();
    void init_PressureStatic();
    void init_PressureDynamic();

    void init_vapour(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &X,
        Array &c, Array &cloud, Array &ice, Array &cloudiness);

    void init_vapour_cloud_ice(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet,
        Array &c, Array &cloud, Array &ice, Array &cloudiness);

    void init_h2s(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, Array &c);

    void OneCategoryIceScheme();

    void Latent_Heat();
    void Forces();

    void steadyQuery();
    void restoreVar(double coeff);



    std::vector<int> im_tropopause; // keep the tropopause layer index
    std::vector<float> m_layer_heights;
    std::vector<double> cloud_loc; // lateral cloudwater distribution
    std::vector<double> r_max; // lateral r_max distribution
    std::vector<double> r_max_add; // lateral r_max distribution
    std::vector<double> t_add; // lateral r_max distribution

    Array_1D rad;
    Array_1D the;
    Array_1D phi;

    Array_2D Topography; // topography
    Array_2D LatentHeat;        // areas of higher latent heat
    Array_2D Precipitation;        // areas of higher precipitation
    Array_2D precipitable_water;// areas of precipitable water in the air
    Array_2D nh3_total;            // areas of higher nh3 concentration
    Array_2D nh3_cloud_total;    // areas of higher nh3_cloud concentration
    Array_2D nh3_ice_total;        // areas of higher nh3_ice concentration
    Array_2D aux_2D_v;            // auxilliar field v
    Array_2D aux_2D_w;            // auxilliar field w
    Array_2D tropopause_height; // local height of the tropopause

    Array t;                    // temperature
    Array u;                    // u-component velocity component in r-direction
    Array v;                    // v-component velocity component in theta-direction
    Array w;                    // w-component velocity component in phi-direction

    Array ch4;                    // methane vapour
    Array ch4_cloud;                // methane cloud
    Array ch4_ice;                    // methane ice
    Array h2o;                    // water vapour
    Array h2o_cloud;                // cloud water
    Array h2o_ice;                    // cloud ice
    Array h2s;                    // water vapour
    Array h2s_cloud;                // cloud water
    Array h2s_ice;                    // cloud ice
    Array nh3;                    // nh3-vapour
    Array nh3_cloud;            // nh3-cloud
    Array nh3_ice;                // nh3-ice
    Array nh4sh;                    // nh4sh-vapour

    Array tn;                    // temperature new
    Array un;                    // u-velocity component in r-direction new
    Array vn;                    // v-velocity component in theta-direction new
    Array wn;                    // w-velocity component in phi-direction new
    Array ch4n;                    // water vapour new
    Array ch4_cloudn;                    // water vapour new
    Array ch4_icen;                // cloud water new
    Array h2on;                    // water vapour new
    Array h2o_cloudn;                    // water vapour new
    Array h2o_icen;                // cloud water new
    Array h2sn;                    // water vapour new
    Array h2s_cloudn;                    // water vapour new
    Array h2s_icen;                // cloud water new
    Array nh3n;                    // nh3 new
    Array nh3_cloudn;            // nh3_cloud new
    Array nh3_icen;                 // nh3_ice new
    Array nh4shn;                    // nh4sh new

    Array rho_mix;         // local mixture density from ideal gas law

    Array massflux_h2s;   // mass flux h2s
    Array massflux_nh3;   // mass flux nh3
    Array massflux_nh4sh;   // mass flux nh4sh

    Array fluxlim_nh4sh;  // TVD flux-limiter correction for nh4sh advection

    Array difflux_h2s;   // diffusive flux h2s
    Array difflux_nh3;   // diffusive flux nh3
    Array difflux_nh4sh;   // diffusive flux nh4sh

    Array thermalmassflux;   // thermal massflux_h2s

    Array cloudiness_ch4; // cloudiness, N in literature
    Array cloudiness_h2o; // cloudiness, N in literature
    Array cloudiness_h2s; // cloudiness, N in literature
    Array cloudiness_nh3; // cloudiness, N in literature

    Array radiation;            // net thermal radiative flux [W/m2]
    Array epsilon;              // layer emissivity
    Array Q_rad;                // radiative heating rate [W/m3]
    Array p_dyn;                // dynamic pressure
    // Previous-iteration dynamic pressure, the n-copy of p_dyn. It did not exist at all
    // until steadyQuery was revived: the routine's pressure case was commented out, its
    // accumulators with it. Maintained by restoreVar with the other n-copies.
    Array p_dynn;               // dynamic pressure, previous iteration
    Array p_stat;                // static pressure

    Array rhs_t;                // auxilliar field RHS temperature
    Array rhs_u;                // auxilliar field RHS u-velocity component
    Array rhs_v;                // auxilliar field RHS v-velocity component
    Array rhs_w;                // auxilliar field RHS w-velocity component

    Array rhs_ch4;                // auxilliar field RHS water vapour
    Array rhs_ch4_cloud;            // auxilliar field RHS cloud water
    Array rhs_ch4_ice;            // auxilliar field RHS cloud water
    Array rhs_h2o;                // auxilliar field RHS water vapour
    Array rhs_h2o_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2o_ice;                // auxilliar field RHS cloud ice
    Array rhs_h2s;                // auxilliar field RHS water vapour
    Array rhs_h2s_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2s_ice;                // auxilliar field RHS cloud ice
    Array rhs_nh3;                // auxilliar field RHS nh3
    Array rhs_nh3_cloud;        // auxilliar field RHS nh3_cloud
    Array rhs_nh3_ice;            // auxilliar field RHS nh3_ice
    Array rhs_nh4sh;                // auxilliar field RHS nh4sh

    Array aux;                // auxilliar field u-velocity component
    Array aux_u;                // auxilliar field u-velocity component
    Array aux_v;                // auxilliar field v-velocity component
    Array aux_w;                // auxilliar field w-velocity component

    Array Q_Latent;                // latent heat
    Array Q_Sensible;            // sensible heat
    Array CoriolisForce;        // coriolis force
    Array CentrifugalForce;             // centrifugal force
    Array BuoyancyForce;        // buoyancy force, Boussinesque approximation
    Array PresGradForce;// pressure gradient force
    Array SeaMount;             // sea mount contour

    Array w_nh3;                // reaction rate nh3
    Array w_h2s;                // reaction rate h2s
    Array w_nh4sh;                // reaction rate nh4sh

    Array j_nh3;                // ordinary-diffusion mass flux of nh3
    Array j_h2s;                // ordinary-diffusion mass flux of h2s
    Array j_nh4sh;                // ordinary-diffusion mass flux of nh4sh

    Array jT_nh3;                // thermo-diffusion mass flux of nh3
    Array jT_h2s;                // thermo-diffusion mass flux of h2s
    Array jT_nh4sh;                // thermo-diffusion mass flux of nh4sh
};
#endif
