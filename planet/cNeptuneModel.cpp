/*
 * Neptune Atmosphere Circulation Model (NACM) applied to laminar flow
 * program for the computation of neptune-atmospherical circulating flows in a spherical shell
 * modeling of the atmosphere with gases, their cloud and ice formation: H2O, H2S, NH3 and NH4SH
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/

#include "cNeptuneModel.h"
#include "TbNept.h"
#include "CwbNept.h"        // ATNEPT_CWB_DIAG, print-only
#include "ConvectiveAdjustmentNept.h"
#include "RadiationNept.h"
#include "TurbulenceNept.h"
#include "PrecipitationNept.h"
#include "PressureSolver.h"
#include "BC_Nept.h"
#include "ChemistryNept.h"
#include "PressureSolverNept.h"
#include "SaturationAdjustmentNept.h"
#include "VelocityInitializerNept.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

// Dry convective adjustment (ConvectiveAdjustmentNept), the SHARED ConvectiveAdjustment<Planet>
// that ATSAT and ATJUP already run. DEFAULT OFF, so every existing ATNEPT run stays
// bit-identical; ATNEPT_CONV_ADJ=1 switches it on. It restores any superadiabatic column to the
// dry adiabat while conserving the column's mass-weighted enthalpy, and nothing else in this
// model does that. Whether Neptune develops such columns at all is unmeasured — switching it on
// and reading the per-iteration report (columns touched, layers mixed, max dT, enthalpy drift)
// is how to find out. Same convention as every other ported module: gated, off, measured later.
// Which pressure solver runs. DEFAULT 0 = ATNEPT's own PressureSolverNept, so every existing
// run stays byte-identical; ATNEPT_PRESS_SOLVER=1 selects the SHARED PressureSolver<Planet> that
// ATSAT and ATJUP run.
//
// THESE ARE NOT THE SAME ALGORITHM, which is why this is a knob rather than a replacement.
// Measured line-for-line after normalising planet names, ATNEPT's solver overlaps the shared one
// by 34.3 % and is half its size (105 significant lines against 218). The shared version carries
// the red-black ordering that makes a sweep mean one thing on any thread count, the obstacle and
// rigid-lid handling, and the metric-radius and coordinate-stretching hooks. What ATNEPT's does
// that the shared one may not is not yet established — that comparison is the next question, and
// it is a physics comparison, not a refactor.
static int press_solver_shared(){
    static const int v = [](){ const char* e = getenv("ATNEPT_PRESS_SOLVER"); return e ? atoi(e) : 0; }();
    return v;
}

// Grey multi-layer radiation (RadiationNept), the SHARED Radiation<Planet> that ATSAT and ATJUP
// already run. DEFAULT ON since 2026-10-11 on the user's word; ATNEPT_RADIATION=0 restores the run
// without it (not byte-checked). It fills the DIAGNOSTIC arrays radiation / epsilon / Q_rad and touches neither t
// nor any rhs — wiring the heating into the temperature equation is a separate step, as it was on
// the other two models.
//
// Neptune is the hard case for this scheme. It absorbs ~1.07 W/m2 of sunlight (S=1.505, A=0.290)
// against an internal flux of 0.433 W/m2, so nearly a third of its budget comes from below —
// where Jupiter, which this scheme was calibrated on, is nearer half but at forty times the
// absolute flux. Whether a grey scheme tuned there behaves at Neptune's temperatures is exactly
// what switching this on is for, and is not established by adding it.
// Turbulence closure (TurbulenceNept), the SHARED Turbulence<Planet> that ATSAT and ATJUP run.
// DEFAULT OFF; ATNEPT_TURB=1 switches it on, and turb_model = "none" also disables it, so both
// have to allow it — on ATSAT those were once independent switches and only one of them decided
// anything, which is worth not repeating. turb_model is k_omega_SST here, as it is on ATSAT.
//
// STAGE ONE OF THREE, following the staging ATSAT used. This fills nue* and the closure
// diagnostics from the velocity field. It does NOT integrate k* and dis* — RungeKutta_Nept
// predates the closure and has no stages for them — and nue* reaches no momentum or scalar
// equation. With the knob ON it moves the diagnostic arrays and nothing else; with it OFF it does
// not run.
//
// ATSAT's warning applies here and has NOT been checked on Neptune: ATSAT's `re` is 1000, which
// made its eddy viscosity ~77x SMALLER than the molecular background — the opposite of ATJUP's
// situation and the opposite of what ATJUP's comment claimed. The first thing to do after
// switching this on is compare nue* against 1/re.
// Precipitation microphysics (PrecipitationNept), the SHARED Precipitation<Planet>. DEFAULT OFF.
// ATSAT runs its precipitation ON by an explicit decision taken after measuring it; ATNEPT's has
// not been measured, and a module that feeds back into the condensate fields is not something to
// switch on by inheritance. It DOES feed back — the condensate it converts is removed from
// cloud/ice in place — so it must run AFTER the SaturationAdjustment calls or they undo it.
// DEFAULT 1 SINCE 2026-10-12 on the user's word ("on Uranus exists no rain of any species"; the
// same on Neptune): the scheme runs, with the rate scale calibrated to the emitted flux
// (Planet::precip_rate_scale()), the corrected t_00 values of NH3, CH4 and H2S (cNeptuneModel.h)
// and the coupling to the species and temperature equations (RHS_Nept_Turb.cpp,
// ATNEPT_PRECIP_COUPLING). The notes above describe the scheme as it stood before that.
// ATNEPT_PRECIP=0 with ATNEPT_T00_LEGACY=1 and ATNEPT_LV_H2S=3.534e6 is the model before (byte-identical
// except for the all-zero Precip_* block in the horizontal VTK).
static int precip_enabled(){
    static const int v = [](){ const char* e = getenv("ATNEPT_PRECIP"); return e ? atoi(e) : 1; }();
    return v;
}

// ATNEPT_TURB: -1 = not set. turb_model IS the switch since 2026-10-11 (the user's word), as on ATJUP
// since 2026-10-09: the closure runs unless turb_model (or ATNEPT_TURB_MODEL, which overrides it) is
// "none" / "laminar" / empty, and param.py's default is k_omega_SST, so a stock run has it ON.
// Until then the closure needed ATNEPT_TURB=1 as well and the configured model alone did nothing.
// ATNEPT_TURB=0 still switches it off; any other value leaves the decision to turb_model.
static int turb_env_enabled(){
    static const int v = [](){ const char* e = getenv("ATNEPT_TURB"); return e ? atoi(e) : -1; }();
    return v;
}

static int radiation_enabled(){
    static const int v = [](){ const char* e = getenv("ATNEPT_RADIATION"); return e ? atoi(e) : 1; }();
    return v;
}

static int conv_adj_enabled(){
    static const int v = [](){ const char* e = getenv("ATNEPT_CONV_ADJ"); return e ? atoi(e) : 0; }();
    return v;
}

cNeptuneModel* cNeptuneModel::m_model = NULL;

const double cNeptuneModel::pi180 = 180.0/M_PI;      // pi180 = 57.3

const double cNeptuneModel::the_degree = 1.0;         // compares to 1° step size laterally
const double cNeptuneModel::phi_degree = 1.0;         // compares to 1° step size longitudinally

const double cNeptuneModel::the0 = 0.0;             // North Pole
const double cNeptuneModel::phi0 = 0.0;             // zero meridian in Greenwich

const double cNeptuneModel::r0 = 1.0; // value much too small, Neptune radius 72000km

const double cNeptuneModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
const double cNeptuneModel::dthe = the_degree/pi180; 
const double cNeptuneModel::dphi = phi_degree/pi180;


cNeptuneModel::cNeptuneModel():
    j_ellipse(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    has_welcome_msg_printed(false){
//    if(PythonStream::is_enable()){
//        backup = std::cout.rdbuf();
//        std::cout.rdbuf(&ps);
//    }
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
    m_model = this;
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
    init_layer_heights();
}

cNeptuneModel::~cNeptuneModel(){
    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }
    m_model = NULL;
}

#include "cNeptuneDefaults.cpp.inc"
/*
*
*/
void cNeptuneModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    if(err){
        doc.PrintError();
        throw std::invalid_argument("   couldn't load config file inside cNeptuneModel");
    }
    XMLElement *atnept = doc.FirstChildElement("atnept");
    if(!atnept){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atnept")->FirstChildElement("common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_neptune = doc.FirstChildElement("atnept")->FirstChildElement("neptune");
    if(!elem_neptune){
        return;
    }
#include "NeptuneLoadConfig.cpp.inc"
}
/*
*
*/
void cNeptuneModel::Run(){
    // THE turbulence gate, resolved once. ATNEPT_TURB and turb_model must BOTH allow it: on ATSAT
    // these were two independent switches and only one of them decided anything, so a run with
    // turb_model = "none" still ran the closure. Both are consulted here and the answer is a bool
    // the rest of the run reads. ATNEPT_TURB_MODEL overrides the configured string first, so the
    // gate is fixed only after the configuration and the override are both final.
    if(const char* tm = getenv("ATNEPT_TURB_MODEL")) turb_model = tm;
    turb_active = (turb_env_enabled() != 0) && !(turb_model.empty() || turb_model == "none" || turb_model == "laminar");


    #ifdef _OPENMP
        printf("\n\n   number of processors: %d\n\n", omp_get_num_procs());

    #pragma omp parallel
        {
            printf("   thread %d of %d in \"Desktop-Dell-XPS 8960\"\n", 
                omp_get_thread_num(), omp_get_num_threads());
        }
    #else
        printf("   OpenMP is not supported\n");
    #endif
        printf("   ended\n\n");

    mkdir(output_path.c_str(), 0777);

    m_model = this;

    cout.precision(6);
    cout.setf(ios::fixed);

    if(!has_welcome_msg_printed)
        print_welcome_msg();

    initMsg();

    panorama_cnt = 0;
    panorama     = 1;
    panorama_step = panorama_print;

    resetArrays();

    dt = 0.001;                                                         //  no dimension

    init_layer_heights();
    TropopauseLocation();
    BC_Nept(*this).initTropopauseLayers();
    VelocityInitializerNept(*this).compute();

    AtomUtils::damp_wiggles(u, nullptr, true, true, true);
    AtomUtils::damp_wiggles(v, nullptr, true, true, true);
    AtomUtils::damp_wiggles(w, nullptr, true, true, true);


    ChemistryNept(*this).ThermalPropertiesNept();
    gam = g * 1.0e3 / cp_mix;                                           // dry adiabatic lapse rate [K/km]: g[m/s²]*1000 / cp[J/(kg·K)]

    // ATNEPT_T_INIT_ANCHOR=1 (cNeptuneModel.h has the note): the reference densities follow the
    // anchored bottom state, one factor on all of them, so that r_mix * R_mix * T_bottom is the
    // bottom pressure. cp_mix, R_mix, mue_mix, k_mix and gam are ratios and stay.
    if(t_init_anchor() != 0){
        const double T_bottom = t_bottom_anchor_K();
        const double f = p_bottom_bar() * 1.0e5 / (R_mix * T_bottom) / r_mix;
        r_h2 *= f; r_he *= f; r_ch4 *= f; r_nh3 *= f; r_h2s *= f; r_h2o *= f; r_nh4sh *= f;
        c_h2 *= f; c_he *= f; c_ch4 *= f; c_nh3 *= f; c_h2s *= f; c_h2o *= f; c_nh4sh *= f;
        r_mix *= f; c_mix *= f;
        const double n = g / (gam * R_ref);
        printf("      ATNEPT: T_INIT_ANCHOR: adiabat through %.1f K at 1 bar, n = %.4f; bottom %.2f K at %.2f bar;"
               " reference densities x %.4f (r_mix %.4f kg/m3); L_atm %.1f km\n",
               t_1bar_K(), n, T_bottom, p_bottom_bar(), f, r_mix, L_atm);
        for(int e = 0; e < 2; e++){
            const double T_tr = e == 0 ? t_equator : t_pole;
            const double h_tr = (T_bottom - T_tr) / gam;
            const double p_tr = p_bottom_bar() * std::pow(T_tr / T_bottom, n);
            const double p_lid = L_atm > h_tr ? p_tr * std::exp(-(L_atm - h_tr) * g / (R_ref * T_tr))
                                              : p_bottom_bar() * std::pow((T_bottom - gam * L_atm) / T_bottom, n);
            printf("      ATNEPT: T_INIT_ANCHOR: %s: isothermal at %.1f K above %.1f km (%.4f bar); lid at %.4e bar\n",
                   e == 0 ? "equator" : "pole   ", T_tr, h_tr, p_tr, p_lid);
        }
    }

    init_temperature();
    init_PressureStatic();

//    goto Printout;

    init_vapour_cloud_ice("CH4", ch4_tropopause, coeff_ch4_A, coeff_ch4_B, 
        coeff_ch4_A_i, coeff_ch4_B_i, t_0_ch4, t_00_ch4, 
        ep_ch4, r_ch4, m_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4,
        ch4, ch4_cloud, ch4_ice, cloudiness_ch4);

    AtomUtils::damp_wiggles(ch4,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("H2O", h2o_tropopause, coeff_h2o_A, coeff_h2o_B, 
        coeff_h2o_A_i, coeff_h2o_B_i, t_0_h2o, t_00_h2o, 
        ep_h2o, r_h2o, m_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o,
        h2o, h2o_cloud, h2o_ice, cloudiness_h2o);

    AtomUtils::damp_wiggles(h2o,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("H2S", h2s_tropopause, coeff_h2s_A, coeff_h2s_B, 
        coeff_h2s_A_i, coeff_h2s_B_i, t_0_h2s, t_00_h2s, 
        ep_h2s, r_h2s, m_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s,
        h2s, h2s_cloud, h2s_ice, cloudiness_h2s);

    AtomUtils::damp_wiggles(h2s,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("NH3", nh3_tropopause, coeff_nh3_A, coeff_nh3_B, 
        coeff_nh3_A_i, coeff_nh3_B_i, t_0_nh3, t_00_nh3, 
        ep_nh3, r_nh3, m_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3,
        nh3, nh3_cloud, nh3_ice, cloudiness_nh3);

    AtomUtils::damp_wiggles(nh3,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_ice,   nullptr, true, true, true);

//    goto Printout;

    SaturationAdjustmentNept(*this).run("H2O",
        coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
        t_0_h2o, t_00_h2o,
        ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
        h2o, h2o_cloud, h2o_ice);

    AtomUtils::damp_wiggles(h2o,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_ice,   nullptr, true, true, true);


    SaturationAdjustmentNept(*this).run("NH3",
        coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
        nh3, nh3_cloud, nh3_ice);

    AtomUtils::damp_wiggles(nh3,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_ice,   nullptr, true, true, true);


    SaturationAdjustmentNept(*this).run("H2S",
        coeff_h2s_A, coeff_h2s_B, coeff_h2s_A_i, coeff_h2s_B_i,
        t_0_h2s, t_00_h2s,
        ep_h2s, lv_h2s, ls_h2s, cp_h2s, r_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, m_h2s,
        h2s, h2s_cloud, h2s_ice);

    AtomUtils::damp_wiggles(h2s,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_ice,   nullptr, true, true, true);


    SaturationAdjustmentNept(*this).run("CH4",
        coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
        ch4, ch4_cloud, ch4_ice);

    AtomUtils::damp_wiggles(ch4,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_ice,   nullptr, true, true, true);


    ChemistryNept(*this).ChemMassRateNept();
    ChemistryNept(*this).FluxLimiterNH4SH();

    AtomUtils::damp_wiggles(massflux_h2s,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(massflux_nh3,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(massflux_nh4sh, nullptr, true, true, true);
    AtomUtils::damp_wiggles(fluxlim_nh4sh,  nullptr, true, true, true);

    ChemistryNept(*this).DiffMassFluxNept();

    AtomUtils::damp_wiggles(difflux_h2s,     nullptr, true, true, true);
    AtomUtils::damp_wiggles(difflux_nh3,     nullptr, true, true, true);
    AtomUtils::damp_wiggles(difflux_nh4sh,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(thermalmassflux, nullptr, true, true, true);

    Forces();
    Latent_Heat();

    init_PressureDynamic();

    BC_Nept(*this).bcTheta();
    BC_Nept(*this).bcPhi();
    BC_Nept(*this).bcRadius();

    restoreVar(1.0);

    // The closure SEEDS k* and dis* algebraically (Turbulence<Planet>::init_fields) before the
    // transport equations start carrying them. It belongs exactly here — after Forces,
    // Latent_Heat, the three boundary passes and the first restoreVar — because init_fields reads
    // the velocity field and the layer geometry, and none of that exists at the top of Run().
    // Putting the call there segfaults; this is ATSAT's slot (cSaturnModel.cpp, right after its
    // own pre-loop restoreVar) and it is the same slot for the same reason.
    if(turb_active) TurbulenceNept(*this).init();

//    goto Printout;



    // ATNEPT_TB_DIAG=<n> (print-only, see TbNept.h): the level-mean temperature after every stage of
    // the loop that can write t, charged to that stage, in K; printed every n iterations.
    static const int tb_every = [](){ const char* e = getenv("ATNEPT_TB_DIAG"); return e ? atoi(e) : 0; }();
    const bool tb_on = tb_every > 0 && im <= NeptTb::MAXLEV;
    NeptTb::st().on = tb_on;
    enum { TB_PRE = 0, TB_SATADJ, TB_PRECIP, TB_PHYS, TB_RK_RESET, TB_RK_REST, TB_BC, TB_CLAMP, TB_AFTER, TB_N };
    static const char* tb_name[TB_N] = { "pressure", "sat.adjust", "precip", "chem/latent", "RK reset", "RK rest",
        "BC", "clamp", "turb/conv" };
    static const char* tb_term[NeptTb::NTERM] = { "RK press.work", "RK transp r", "RK transp the", "RK transp phi",
        "RK diffusion", "RK th.massflux", "RK radiation", "RK precip" };
    std::vector<std::vector<double>> tb_sum(TB_N, std::vector<double>(im, 0.0)), tb_tsum(NeptTb::NTERM, std::vector<double>(im, 0.0));
    std::vector<double> tb_last(im, 0.0), tb_first(im, 0.0);
    double tb_norm = 0.0;
    auto tb_lev = [&](Array& a) -> std::vector<double> {
        std::vector<double> m_(im, 0.0);
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < im; i++){
            double s_ = 0.0;
            for(int j = 0; j < jm; j++){
                const double sj = std::sin(the.z[j]);
                for(int k = 0; k < km; k++) s_ += sj * a.x[i][j][k];
            }
            m_[i] = s_ * t_ref / tb_norm;
        }
        return m_; };
    auto tb_mark = [&](int stage){
        if(!tb_on) return;
        const std::vector<double> w_ = tb_lev(t);
        for(int i = 0; i < im; i++){ tb_sum[stage][i] += w_[i] - tb_last[i]; tb_last[i] = w_[i]; } };
    if(tb_on){
        for(int j = 0; j < jm; j++) tb_norm += std::sin(the.z[j]) * km;
        tb_last = tb_first = tb_lev(t);
    }

    // ATNEPT_CWB_DIAG=<n> (print-only, see CwbNept.h): the H2O, NH3, CH4 and H2S columns after
    // every stage of the loop that can write them, charged to that stage, in g/m2; printed every
    // n iterations. The stages are the ones ATNEPT_TB_DIAG marks.
    static const int cwb_every = [](){ const char* e = getenv("ATNEPT_CWB_DIAG"); return e ? atoi(e) : 0; }();
    const bool cwb_on = cwb_every > 0;
    NeptCwb::st().on = cwb_on;
    static const char* cw_gas[NeptCwb::NGAS] = { "H2O", "NH3", "CH4", "H2S" };
    static const char* cw_name[TB_N] = { "before the adjustment (pressure)", "SaturationAdjustment (all calls)",
        "Precipitation (in place)", "scalar sync, chemistry, forces, radiation", "RungeKutta: reset to the n-level copies",
        "RungeKutta: floor, limiters and remainder", "boundary conditions", "clampNegativeSpecies", "turbulence, convective adjustment, filters" };
    static const char* cw_term[NeptCwb::NTERM] = { "RungeKutta: radial transport   -u dq/dr", "RungeKutta: meridional transport",
        "RungeKutta: zonal transport", "RungeKutta: diffusion", "RungeKutta: chemistry source", "RungeKutta: precipitation source",
        "(reference, not a term) q*div(u)" };
    Array* cw_f[NeptCwb::NGAS][3]  = { { &h2o,  &h2o_cloud,  &h2o_ice  }, { &nh3,  &nh3_cloud,  &nh3_ice  }, { &ch4,  &ch4_cloud,  &ch4_ice  }, { &h2s,  &h2s_cloud,  &h2s_ice  } };
    Array* cw_fn[NeptCwb::NGAS][3] = { { &h2on, &h2o_cloudn, &h2o_icen }, { &nh3n, &nh3_cloudn, &nh3_icen }, { &ch4n, &ch4_cloudn, &ch4_icen }, { &h2sn, &h2s_cloudn, &h2s_icen } };
    double cw_sum[NeptCwb::NGAS][TB_N] = {{0.0}}, cw_tsum[NeptCwb::NGAS][NeptCwb::NTERM] = {{0.0}};
    double cw_norm = 0.0, cw_last[NeptCwb::NGAS] = {0.0}, cw_first[NeptCwb::NGAS] = {0.0}, cw_first0[NeptCwb::NGAS] = {0.0};
    // i_hi = 0: every level with a thickness; i_hi = 1: the plane i = 0 alone.
    auto cwb_col = [&](Array* f[3], int i_hi) -> double {
        double s_ = 0.0;
        const int ie = i_hi > 0 ? i_hi : im;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:s_)
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                const double sj = std::sin(the.z[j]);
                for(int i = 0; i < ie; i++){
                    const double dz = layer_thickness_m(i);
                    if(dz > 0.0) s_ += sj * dz * (f[0]->x[i][j][k] + f[1]->x[i][j][k] + f[2]->x[i][j][k]);
                }
            }
        }
        return s_ * 1.0e3 / cw_norm; };
    auto cwb_mark = [&](int stage){
        if(!cwb_on) return;
        for(int g_ = 0; g_ < NeptCwb::NGAS; g_++){
            const double w_ = cwb_col(cw_f[g_], 0);
            cw_sum[g_][stage] += w_ - cw_last[g_]; cw_last[g_] = w_;
        } };
    if(cwb_on){
        for(int j = 0; j < jm; j++) cw_norm += std::sin(the.z[j]) * km;
        for(int g_ = 0; g_ < NeptCwb::NGAS; g_++){
            cw_last[g_] = cw_first[g_] = cwb_col(cw_f[g_], 0);
            cw_first0[g_] = cwb_col(cw_f[g_], 1);
        }
    }

    for(iter_n = 1; iter_n <= nm; iter_n++){

        auto begin = std::chrono::high_resolution_clock::now();

        cout << endl << endl;
        cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        cout << " 3D Neptune iterational process" << endl;
        cout << " present state of the computation " << endl << endl
             << " ======================== iteration number n = " << iter_n << " ========================= " << endl << endl
             << "    max total iteration number nm = " << nm << endl
             << "    checkpoint when to write 3D-panorama = " << checkpoint << endl
             << "    panorama_print = " << panorama_print << endl << endl;

        if(iter_n % 2 == 0){

            if(press_solver_shared()) PressureSolver<cNeptuneModel>(*this).run();
            else                      PressureSolverNept(*this).run();
            AtomUtils::damp_wiggles(p_dyn, nullptr, true, true, true);
            tb_mark(TB_PRE); cwb_mark(TB_PRE);

            SaturationAdjustmentNept(*this).run("H2O",
                coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
                t_0_h2o, t_00_h2o,
                ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
                C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
                h2o, h2o_cloud, h2o_ice);

            SaturationAdjustmentNept(*this).run("H2S",
                coeff_h2s_A, coeff_h2s_B, coeff_h2s_A_i, coeff_h2s_B_i,
                t_0_h2s, t_00_h2s,
                ep_h2s, lv_h2s, ls_h2s, cp_h2s, r_h2s,
                C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, m_h2s,
                h2s, h2s_cloud, h2s_ice);

            SaturationAdjustmentNept(*this).run("NH3",
                coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
                t_0_nh3, t_00_nh3,
                ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
                C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
                nh3, nh3_cloud, nh3_ice);

            SaturationAdjustmentNept(*this).run("CH4",
                coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
                t_0_ch4, t_00_ch4,
                ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
                C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
                ch4, ch4_cloud, ch4_ice);

            // Must follow the saturation adjustments: it removes condensate in place.
            tb_mark(TB_SATADJ); cwb_mark(TB_SATADJ);
            if(precip_enabled()) PrecipitationNept(*this).run();
            tb_mark(TB_PRECIP); cwb_mark(TB_PRECIP);

            // ATNEPT_RK_SCALAR_SYNC=<0|1|2> (2026-10-11; DEFAULT 2 on the user's word, 0 restores):
            // ATJUP's repair of 2026-10-08, as ported to ATSAT the same day. RungeKuttaNept
            // integrates every scalar from its n-level copy (h2on, ..., tn), and restoreVar
            // refreshed those at the END of the previous iteration -- before the saturation
            // adjustment above ran. So everything the adjustment writes directly (the condensate
            // it forms, the latent heat it releases) was overwritten by the step that follows.
            // 1 copies the twelve condensable fields into their n-level copies, so the step starts
            // from the adjusted state; 2 copies the temperature too, so the latent heat goes with
            // the condensate. The in-place precipitation above is kept by it as well.
            static const int rk_scalar_sync = [](){ const char* e = getenv("ATNEPT_RK_SCALAR_SYNC"); return e ? atoi(e) : 2; }();
            if(rk_scalar_sync >= 1){
                #pragma omp parallel for collapse(2) schedule(static)
                for(int i = 0; i < im; i++){
                    for(int j = 0; j < jm; j++){
                        for(int k = 0; k < km; k++){
                            h2on.x[i][j][k]       = h2o.x[i][j][k];
                            h2o_cloudn.x[i][j][k] = h2o_cloud.x[i][j][k];
                            h2o_icen.x[i][j][k]   = h2o_ice.x[i][j][k];
                            h2sn.x[i][j][k]       = h2s.x[i][j][k];
                            h2s_cloudn.x[i][j][k] = h2s_cloud.x[i][j][k];
                            h2s_icen.x[i][j][k]   = h2s_ice.x[i][j][k];
                            nh3n.x[i][j][k]       = nh3.x[i][j][k];
                            nh3_cloudn.x[i][j][k] = nh3_cloud.x[i][j][k];
                            nh3_icen.x[i][j][k]   = nh3_ice.x[i][j][k];
                            ch4n.x[i][j][k]       = ch4.x[i][j][k];
                            ch4_cloudn.x[i][j][k] = ch4_cloud.x[i][j][k];
                            ch4_icen.x[i][j][k]   = ch4_ice.x[i][j][k];
                            if(rk_scalar_sync >= 2) tn.x[i][j][k] = t.x[i][j][k];
                        }
                    }
                }
            }

            ChemistryNept(*this).DiffMassFluxNept();

            AtomUtils::damp_wiggles(difflux_h2s,     nullptr, true, true, true);
            AtomUtils::damp_wiggles(difflux_nh3,     nullptr, true, true, true);
            AtomUtils::damp_wiggles(difflux_nh4sh,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(thermalmassflux, nullptr, true, true, true);

            ChemistryNept(*this).ChemMassRateNept();
            ChemistryNept(*this).FluxLimiterNH4SH();

            AtomUtils::damp_wiggles(massflux_h2s,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(massflux_nh3,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(massflux_nh4sh, nullptr, true, true, true);
            AtomUtils::damp_wiggles(fluxlim_nh4sh,  nullptr, true, true, true);

            Forces();
            Latent_Heat();
            tb_mark(TB_PHYS); cwb_mark(TB_PHYS);

        }  // if loop

        std::vector<double> tb_n;
        if(tb_on){
            tb_n = tb_lev(tn);                                          // what the step starts from
            for(int t_ = 0; t_ < NeptTb::MAXTHR; t_++) for(int q_ = 0; q_ < NeptTb::NTERM; q_++)
                for(int i = 0; i < im; i++) NeptTb::st().acc[t_][q_][i] = 0.0;
        }

        double cw_n[NeptCwb::NGAS] = {0.0};
        if(cwb_on){
            for(int g_ = 0; g_ < NeptCwb::NGAS; g_++) cw_n[g_] = cwb_col(cw_fn[g_], 0);    // what the step starts from
            for(int t_ = 0; t_ < NeptCwb::MAXTHR; t_++) for(int g_ = 0; g_ < NeptCwb::NGAS; g_++)
                for(int q_ = 0; q_ < 8; q_++) NeptCwb::st().acc[t_][g_][q_] = 0.0;
        }

        RungeKuttaNept();

        if(cwb_on){
            for(int g_ = 0; g_ < NeptCwb::NGAS; g_++){
                const double w_ = cwb_col(cw_f[g_], 0);
                double step_terms = 0.0;
                for(int q_ = 0; q_ < NeptCwb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < NeptCwb::MAXTHR; t_++) a_ += NeptCwb::st().acc[t_][g_][q_];
                    a_ *= dt / 6.0 * 1.0e3 / cw_norm;
                    cw_tsum[g_][q_] += a_;
                    if(q_ != NeptCwb::T_QDIV) step_terms += a_;
                }
                cw_sum[g_][TB_RK_RESET] += cw_n[g_] - cw_last[g_];
                cw_sum[g_][TB_RK_REST]  += w_ - cw_n[g_] - step_terms;
                cw_last[g_] = w_;
            }
        }

        if(tb_on){
            const std::vector<double> w_ = tb_lev(t);
            for(int i = 0; i < im; i++){
                double step_terms = 0.0;
                for(int q_ = 0; q_ < NeptTb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < NeptTb::MAXTHR; t_++) a_ += NeptTb::st().acc[t_][q_][i];
                    a_ *= dt / 6.0 * t_ref / tb_norm;
                    tb_tsum[q_][i] += a_; step_terms += a_;
                }
                tb_sum[TB_RK_RESET][i] += tb_n[i] - tb_last[i];
                tb_sum[TB_RK_REST][i]  += w_[i] - tb_n[i] - step_terms;
                tb_last[i] = w_[i];
            }
        }

        // ATNEPT_DAMP_T_VERT (2026-10-11, mirrors ATURAN_DAMP_T_VERT): the strength of the VERTICAL
        // pass of the Shapiro 1-2-1 filter on t; the horizontal passes are untouched. damp_wiggles
        // applies its axes in the order k, j, i, so the split below keeps the order of the single
        // (true, true, true) call it replaces. 1.0 is that call; 0.0 drops the vertical pass.
        // WHY: ATNEPT_TB_DIAG, 224 iterations (satchk/giants/NEPTtbA), level-mean change at levels
        // 1 / 5 / 20 / 35 / 39: -132.1 / -77.1 / +2.3 / +62.3 / +127.7 K, of which this stage
        // (the filter and the boundary conditions) -130.9 / -76.7 / +0.7 / +96.4 / +127.2 and
        // every term of rhs_t together under 4 K except the thermal mass flux at level 35 (-35).
        // The filter's ends are clamped, so on a straight lapse it pulls both ends towards their
        // neighbours every iteration and the flattening works its way inwards: a noise filter
        // acting as the model's vertical heat transport. ATJUP and ATSAT do not filter t at all.
        // DEFAULT 0.0 SINCE 2026-10-11 on the user's word (1.0 restores the filter of every run
        // before that). MEASURED with ATNEPT_TB_DIAG, 224 iterations, level-mean change in K at levels
        // 1 / 20 / 35 / 39 (satchk/giants/NEPTtbB against NEPTtbE):
        //     vertical pass 1.0, thermal mass flux 1     -132.1 / +2.3 / +62.3 / +127.7
        //     vertical pass 0.0, thermal mass flux 0     -5.5 / -1.1 / +1.1 / +7.6
        // mean OLR 94.9 -> 5.98 W/m2 against 0.700 in, photosphere 203.1 -> 83.3 K; exit 0, no non-finite
        // value. With the pass off and the thermal mass flux still on, that term digs a cold hole
        // (-51 K at level 35), which is why ATNEPT_THERMAL_MASSFLUX went to 0 in the same change.
        // What is left is the diffusion term at the two radial ends (the 2-point wall value does
        // not continue a straight lapse) and a rougher profile: the README's sweep measured the
        // mean second difference 3.1 times larger at 0.0 than at 1.0, with nothing oscillating.
        static const double damp_t_vert = [](){
            const char* e = getenv("ATNEPT_DAMP_T_VERT"); return e ? atof(e) : 0.0; }();
        AtomUtils::damp_wiggles(t, nullptr, false, true, true);              // horizontal, unchanged
        if(damp_t_vert != 0.0)
            AtomUtils::damp_wiggles(t, nullptr, true, false, false, damp_t_vert);
        AtomUtils::damp_wiggles(u, nullptr, true, true, true);
        AtomUtils::damp_wiggles(v, nullptr, true, true, true);
        AtomUtils::damp_wiggles(w, nullptr, true, true, true);

        BC_Nept(*this).bcRadius();                                      // extrapolation in i-direction along grid boundaries
        BC_Nept(*this).bcTheta();                                       // extrapolation in j-direction along grid boundaries
        BC_Nept(*this).bcPhi();                                         // extrapolation in k-direction along grid boundaries
        tb_mark(TB_BC); cwb_mark(TB_BC);
        clampNegativeSpecies();                                         // ATNEPT_SPECIES_CLAMP, default on (BC_Nept.h)
        tb_mark(TB_CLAMP); cwb_mark(TB_CLAMP);

        // How far the run is from a steady state, and WHERE. MUST run BEFORE restoreVar: it
        // differences each field against the n-copy restoreVar is about to overwrite, so after
        // that call every number it prints is identically zero. That is the reason it was never
        // called from anywhere useful, and the reason it is called here.
        //
        // Same cadence as printMinMax. ATNEPT_STEADY=0 switches it off. This ADDS LINES TO THE
        // LOG that ATNEPT has never printed — the routine has never run in this model — but it
        // writes nothing and changes no output file.
        static const int steady_on = [](){
            const char* e = getenv("ATNEPT_STEADY"); return e ? atoi(e) : 1; }();
        if(steady_on && iter_n % checkpoint == 0) steadyQuery();

        restoreVar(1.0);

        // After the state has been advanced and the boundaries applied: put any
        // superadiabatic column back on the dry adiabat. Off by default (ATNEPT_CONV_ADJ).
        if(turb_active) TurbulenceNept(*this).run();
        if(radiation_enabled()) RadiationNept(*this).run();
        if(conv_adj_enabled()) ConvectiveAdjustmentNept(*this).run();
        tb_mark(TB_AFTER); cwb_mark(TB_AFTER);
        if(tb_on && iter_n % tb_every == 0){
            static const int lv_[] = { 0, 1, 2, 3, 5, 8, 12, 16, 20, 24, 28, 32, 35, 37, 38, 39, 40 };
            const int nl_ = (int)(sizeof(lv_) / sizeof(lv_[0]));
            printf("\n      ATNEPT: [TB] level-mean temperature, cumulative over iterations %d..%d, K\n      ATNEPT: [TB] %-15s", 1, iter_n, "level");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %8d", lv_[l]);
            printf("\n      ATNEPT: [TB] %-15s", "T now");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %8.2f", tb_last[lv_[l]]);
            printf("\n      ATNEPT: [TB] %-15s", "change");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_last[lv_[l]] - tb_first[lv_[l]]);
            for(int q_ = 0; q_ < TB_N; q_++){
                printf("\n      ATNEPT: [TB] %-15s", tb_name[q_]);
                for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_sum[q_][lv_[l]]);
                if(q_ == TB_RK_RESET) for(int r_ = 0; r_ < NeptTb::NTERM; r_++){
                    printf("\n      ATNEPT: [TB] %-15s", tb_term[r_]);
                    for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_tsum[r_][lv_[l]]);
                }
            }
            printf("\n");
        }

        if(cwb_on && iter_n % cwb_every == 0){
            for(int g_ = 0; g_ < NeptCwb::NGAS; g_++){
                double sum_ = 0.0;
                printf("\n      ATNEPT: [CWB] %s column (vapour+cloud+ice), cumulative over iterations 1..%d, g/m2:  start %.9e  now %.9e  change %+.6e (%+.4f %%)\n",
                       cw_gas[g_], iter_n, cw_first[g_], cw_last[g_], cw_last[g_] - cw_first[g_], 100.0 * (cw_last[g_] / cw_first[g_] - 1.0));
                const double p0_ = cwb_col(cw_f[g_], 1) - cw_first0[g_], in_ = (cw_last[g_] - cw_first[g_]) - p0_;
                printf("      ATNEPT: [CWB] %s   of the change: plane i=0 (set by bcRadius) %+.6e, integrated cells i>=1 %+.6e (%+.4f %% of their start)\n",
                       cw_gas[g_], p0_, in_, 100.0 * in_ / (cw_first[g_] - cw_first0[g_]));
                for(int q_ = 0; q_ < TB_N; q_++){
                    printf("      ATNEPT: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], cw_name[q_], cw_sum[g_][q_]); sum_ += cw_sum[g_][q_];
                    if(q_ == TB_RK_RESET)
                        for(int r_ = 0; r_ < NeptCwb::NTERM; r_++){
                            printf("      ATNEPT: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], cw_term[r_], cw_tsum[g_][r_]);
                            if(r_ != NeptCwb::T_QDIV) sum_ += cw_tsum[g_][r_];
                        }
                }
                printf("      ATNEPT: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], "unattributed", cw_last[g_] - cw_first[g_] - sum_);
            }
        }

        panorama_cnt++;

        if(iter_n % checkpoint == 0){
            printMinMax();
            writeData();
        }

        if(panorama_cnt == panorama_print) panorama_cnt = 0;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for one time step\n", elapsed.count() * 1e-9);

    }  // end for iter_n

    clampNegativeReport();

    cout << endl << "      ATNEPT: run_3D_loop ended ..........................." << endl;

/*
Printout:
    printMinMax();
    writeData();
    print_final_msg();
*/

    return;
}
/*
*
*/
void cNeptuneModel::resetArrays(){
    cout << endl << "      ATNEPT: resetArrays" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    Topography.initArray_2D(jm, km, 0.0); // topography
    LatentHeat.initArray_2D(jm, km, 0.0);            // areas of higher latent heat
    Precipitation.initArray_2D(jm, km, 0.0);         // areas of higher precipitation
    precipitable_water.initArray_2D(jm, km, 0.0);    // areas of precipitable water in the air
    nh3_total.initArray_2D(jm, km, 0.0);             // areas of higher nh3 concentration
    nh3_cloud_total.initArray_2D(jm, km, 0.0);       // areas of higher nh3_cloud concentration
    nh3_ice_total.initArray_2D(jm, km, 0.0);         // areas of higher nh3_ice concentration
    aux_2D_v.initArray_2D(jm, km, 0.0);              // auxilliar field v
    aux_2D_w.initArray_2D(jm, km, 0.0);              // auxilliar field w
    tropopause_height.initArray_2D(jm, km, 0.0); // local height of the tropopause

    t.initArray(im, jm, km, ta);                    // temperature
    u.initArray(im, jm, km, ua);                    // u-component velocity component in r-direction
    v.initArray(im, jm, km, va);                    // v-component velocity component in theta-direction
    w.initArray(im, jm, km, wa);                    // w-component velocity component in phi-direction

    h2o.initArray(im, jm, km, 0.0);                    // water vapour
    h2o_cloud.initArray(im, jm, km, 0.0);                // cloud water
    h2o_ice.initArray(im, jm, km, 0.0);                  // cloud ice
    ch4.initArray(im, jm, km, 0.0);                    // methane vapour
    ch4_cloud.initArray(im, jm, km, 0.0);                // methane cloud
    ch4_ice.initArray(im, jm, km, 0.0);                  // methane ice
    h2s.initArray(im, jm, km, 0.0);                    // water vapour
    h2s_cloud.initArray(im, jm, km, 0.0);                // cloud water
    h2s_ice.initArray(im, jm, km, 0.0);                  // cloud ice
    nh3.initArray(im, jm, km, 0.0);                 // nh3-vapour
    nh3_cloud.initArray(im, jm, km, 0.0);           // nh3-cloud
    nh3_ice.initArray(im, jm, km, 0.0);             // nh3-ice
    nh4sh.initArray(im, jm, km, 0.0);                 // nh4sh-vapour

    tn.initArray(im, jm, km, ta);                    // temperature new
    un.initArray(im, jm, km, ua);                    // u-velocity component in r-direction new
    vn.initArray(im, jm, km, va);                    // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, wa);                    // w-velocity component in phi-direction new

    ch4n.initArray(im, jm, km, 0.0);                    // methane vapour new
    ch4_cloudn.initArray(im, jm, km, 0.0);                // methane cloud new
    ch4_icen.initArray(im, jm, km, 0.0);                    // methane ice new
    h2on.initArray(im, jm, km, 0.0);                    // water vapour new
    h2o_cloudn.initArray(im, jm, km, 0.0);                // cloud water new
    h2o_icen.initArray(im, jm, km, 0.0);                    // cloud ice new
    h2sn.initArray(im, jm, km, 0.0);                    // water vapour new
    h2s_cloudn.initArray(im, jm, km, 0.0);                // cloud water new
    h2s_icen.initArray(im, jm, km, 0.0);                    // cloud ice new
    nh3n.initArray(im, jm, km, 0.0);                // nh3 new
    nh3_cloudn.initArray(im, jm, km, 0.0);            // nh3_cloud new
    nh3_icen.initArray(im, jm, km, 0.0);            // nh3_ice new
    nh4shn.initArray(im, jm, km, 0.0);                // nh4sh new

    massflux_h2s.initArray(im, jm, km, 0.0);   // mass flux h2s
    massflux_nh3.initArray(im, jm, km, 0.0);   // mass flux nh3
    massflux_nh4sh.initArray(im, jm, km, 0.0);   // mass flux nh4sh

    fluxlim_nh4sh.initArray(im, jm, km, 0.0);  // TVD flux-limiter correction

    difflux_h2s.initArray(im, jm, km, 0.0);   // diffusive flux h2s
    difflux_nh3.initArray(im, jm, km, 0.0);   // diffusive flux nh3
    difflux_nh4sh.initArray(im, jm, km, 0.0); // diffusive flux nh4sh

    thermalmassflux.initArray(im, jm, km, 0.0);   // thermal massflux_h2s

    acc_tke.initArray(im, jm, km, 0.0);
    acc_dis.initArray(im, jm, km, 0.0);
    rhs_tke.initArray(im, jm, km, 0.0);
    rhs_dis.initArray(im, jm, km, 0.0);
    acc_t.initArray(im, jm, km, 0.0);
    acc_u.initArray(im, jm, km, 0.0);
    acc_v.initArray(im, jm, km, 0.0);
    acc_w.initArray(im, jm, km, 0.0);
    acc_ch4.initArray(im, jm, km, 0.0);
    acc_ch4_cloud.initArray(im, jm, km, 0.0);
    acc_ch4_ice.initArray(im, jm, km, 0.0);
    acc_h2o.initArray(im, jm, km, 0.0);
    acc_h2o_cloud.initArray(im, jm, km, 0.0);
    acc_h2o_ice.initArray(im, jm, km, 0.0);
    acc_h2s.initArray(im, jm, km, 0.0);
    acc_h2s_cloud.initArray(im, jm, km, 0.0);
    acc_h2s_ice.initArray(im, jm, km, 0.0);
    acc_nh3.initArray(im, jm, km, 0.0);
    acc_nh3_cloud.initArray(im, jm, km, 0.0);
    acc_nh3_ice.initArray(im, jm, km, 0.0);
    acc_nh4sh.initArray(im, jm, km, 0.0);
    P_rain.initArray(im, jm, km, 0.0);
    P_snow.initArray(im, jm, km, 0.0);
    P_graupel.initArray(im, jm, km, 0.0);
    P_nh3_rain.initArray(im, jm, km, 0.0);
    P_nh3_snow.initArray(im, jm, km, 0.0);
    P_nh3_graupel.initArray(im, jm, km, 0.0);
    P_ch4_rain.initArray(im, jm, km, 0.0);
    P_ch4_snow.initArray(im, jm, km, 0.0);
    P_ch4_graupel.initArray(im, jm, km, 0.0);
    P_nh4sh.initArray(im, jm, km, 0.0);
    Q_precip.initArray(im, jm, km, 0.0);
    S_precip_h2o.initArray(im, jm, km, 0.0);
    S_precip_h2o_cloud.initArray(im, jm, km, 0.0);
    S_precip_h2o_ice.initArray(im, jm, km, 0.0);
    S_precip_nh3.initArray(im, jm, km, 0.0);
    S_precip_nh3_cloud.initArray(im, jm, km, 0.0);
    S_precip_nh3_ice.initArray(im, jm, km, 0.0);
    S_precip_ch4.initArray(im, jm, km, 0.0);
    S_precip_ch4_cloud.initArray(im, jm, km, 0.0);
    S_precip_ch4_ice.initArray(im, jm, km, 0.0);
    precip_srf_h2o.initArray_2D(jm, km, 0.0);
    precip_srf_nh3.initArray_2D(jm, km, 0.0);
    precip_srf_ch4.initArray_2D(jm, km, 0.0);
    P_h2s_rain.initArray(im, jm, km, 0.0);
    P_h2s_snow.initArray(im, jm, km, 0.0);
    P_h2s_graupel.initArray(im, jm, km, 0.0);
    S_precip_h2s.initArray(im, jm, km, 0.0);
    S_precip_h2s_cloud.initArray(im, jm, km, 0.0);
    S_precip_h2s_ice.initArray(im, jm, km, 0.0);
    precip_srf_h2s.initArray_2D(jm, km, 0.0);
    precip_srf_nh4sh.initArray_2D(jm, km, 0.0);
    precip_srf_total.initArray_2D(jm, km, 0.0);
    tke.initArray(im, jm, km, 0.0);
    dis.initArray(im, jm, km, 0.0);
    tken.initArray(im, jm, km, 0.0);
    disn.initArray(im, jm, km, 0.0);
    nue.initArray(im, jm, km, 0.0);
    nue_t.initArray(im, jm, km, 0.0);
    prod.initArray(im, jm, km, 0.0);
    tke_source.initArray(im, jm, km, 0.0);
    dis_source.initArray(im, jm, km, 0.0);
    vel_star.initArray_2D(jm, km, 0.0);
    radiation.initArray(im, jm, km, 0.0);            // net thermal radiative flux [W/m2]
    epsilon.initArray(im, jm, km, 0.0);              // layer emissivity
    Q_rad.initArray(im, jm, km, 0.0);                // radiative heating rate [W/m3]
    p_dyn.initArray(im, jm, km, pa);                // dynamic pressure
    p_hydro.initArray(im, jm, km, 0.0);             // hydrostatic pressure perturbation
    p_dynn.initArray(im, jm, km, pa);               // dynamic pressure, previous iteration
    p_stat.initArray(im, jm, km, pa);                // static pressure
    rho_mix.initArray(im, jm, km, 0.0);          // local mixture density


    rhs_t.initArray(im, jm, km, 0.0);                // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.0);                // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.0);                // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.0);                // auxilliar field RHS w-velocity component

    rhs_ch4.initArray(im, jm, km, 0.0);                // auxilliar field RHS methane vapour
    rhs_ch4_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS methane cloud
    rhs_ch4_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS methane ice
    rhs_h2o.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_h2o_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS cloud water
    rhs_h2o_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS cloud ice
    rhs_h2s.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_h2s_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS cloud water
    rhs_h2s_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS cloud ice
    rhs_nh3.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh3
    rhs_nh3_cloud.initArray(im, jm, km, 0.0);        // auxilliar field RHS nh3_cloud
    rhs_nh3_ice.initArray(im, jm, km, 0.0);            // auxilliar field RHS nh3_ice
    rhs_nh4sh.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh4sh

    aux.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_u.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.0);                // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.0);                // auxilliar field w-velocity component

    Q_Latent.initArray(im, jm, km, 0.0);                // latent heat
    Q_Sensible.initArray(im, jm, km, 0.0);            // sensible heat
    CoriolisForce.initArray(im, jm, km, 0.0);        // coriolis force
    CentrifugalForce.initArray(im, jm, km, 0.0);             // centrifugal force
    BuoyancyForce.initArray(im, jm, km, 0.0);        // buoyancy force, Boussinesque approximation
    PresGradForce.initArray(im, jm, km, 0.0);// pressure gradient force
    SeaMount.initArray(im, jm, km, 0.0);             // sea mount contour

    cloudiness_ch4.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_h2o.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_h2s.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_nh3.initArray(im, jm, km, 0.0); // cloudiness, N in literature

    w_nh3.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh3
    w_h2s.initArray(im, jm, km, 0.0);                // chemical reaction rate of h2s
    w_nh4sh.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh4sh

    j_nh3.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of nh3
    j_h2s.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of h2s
    j_nh4sh.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of nh4sh

    jT_nh3.initArray(im, jm, km, 0.0);                // thermo-diffusion mass flux of nh3
    jT_h2s.initArray(im, jm, km, 0.0);                // thermo-diffusion mass flux of h2s
    jT_nh4sh.initArray(im, jm, km, 0.0);                // thermo-diffusion mass flux of nh4sh

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for resetArrays\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: resetArrays ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::restoreVar(double coeff){
//    cout << endl << "      ATNEPT: restoreVar" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                // p_dynn is the previous-iteration dynamic pressure. Nothing maintained it
                // and nothing allocated it, while steadyQuery's pressure case sat commented
                // out because of that. It belongs with the other n-copies.
                p_dynn.x[i][j][k] = coeff * p_dyn.x[i][j][k];
                // k* and dis* take the same start-of-step copies the other prognostic fields
                // get; without them the RK4 stages would integrate from a moving base. With the
                // closure off both fields are identically zero, so these are zero too.
                tken.x[i][j][k] = coeff * tke.x[i][j][k];
                disn.x[i][j][k] = coeff * dis.x[i][j][k];
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                un.x[i][j][k] = coeff * u.x[i][j][k];
                vn.x[i][j][k] = coeff * v.x[i][j][k];
                wn.x[i][j][k] = coeff * w.x[i][j][k];
                ch4n.x[i][j][k] = coeff * ch4.x[i][j][k];
                ch4_cloudn.x[i][j][k] = coeff * ch4_cloud.x[i][j][k];
                ch4_icen.x[i][j][k] = coeff * ch4_ice.x[i][j][k];
                h2on.x[i][j][k] = coeff * h2o.x[i][j][k];
                h2o_cloudn.x[i][j][k] = coeff * h2o_cloud.x[i][j][k];
                h2o_icen.x[i][j][k] = coeff * h2o_ice.x[i][j][k];
                h2sn.x[i][j][k] = coeff * h2s.x[i][j][k];
                h2s_cloudn.x[i][j][k] = coeff * h2s_cloud.x[i][j][k];
                h2s_icen.x[i][j][k] = coeff * h2s_ice.x[i][j][k];
                nh3n.x[i][j][k] = coeff * nh3.x[i][j][k];
                nh3_cloudn.x[i][j][k] = coeff * nh3_cloud.x[i][j][k];
                nh3_icen.x[i][j][k] = coeff * nh3_ice.x[i][j][k];
                nh4shn.x[i][j][k] = coeff * nh4sh.x[i][j][k];
//                nh4sh_cloudn.x[i][j][k] = coeff * nh4sh_cloud.x[i][j][k];
            }
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for restoreVar\n", elapsed.count() * 1e-9);

//    cout << "      ATNEPT: restoreVar ended" << endl;

    return;
}
/*
*
*/
/*
    cout << endl << "      ATNEPT: thermodynamic fft_gaussian_filter_3d in Run begin ......................." << endl;

    fft_gaussian_filter_3d(nh3,1);
    fft_gaussian_filter_3d(nh3_cloud,1);
    fft_gaussian_filter_3d(nh3_ice,1);

    fft_gaussian_filter_3d(nh4sh,1);
//    fft_gaussian_filter_3d(nh4sh_cloud,1);
//    fft_gaussian_filter_3d(nh4sh_ice,1);

    fft_gaussian_filter_3d(h2s,1);
    fft_gaussian_filter_3d(h2s_cloud,1);
    fft_gaussian_filter_3d(h2s_ice,1);

    fft_gaussian_filter_3d(h2o,1);
    fft_gaussian_filter_3d(h2o_cloud,1);
    fft_gaussian_filter_3d(h2o_ice,1);

    cout << endl << "      ATNEPT: thermodynamic fft_gaussian_filter_3d in Run end ......................." << endl;
*/


