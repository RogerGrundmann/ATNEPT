/*
 * Atmosphere General Circulation Modell(ATNEPT) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cNeptuneModel.h"
#include "SaturationAdjustmentNept.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

void cNeptuneModel::init_temperature(){
    cout << endl << "      ATNEPT: init_Temperature" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    int j_half = (jm-1)/2;
    double d_j_half = (double)j_half;
    double t_eff = t_pole - t_equator; // non-dimensional
    double d_j = 0.0;
    double height = 0.0;

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            d_j = (double)(j);
            // t_equator and t_pole are the TOP (i=im-1) temperatures [K].
            // T_bottom = T_top + gam * L_atm; gam [K/km], L_atm [km] via get_layer_height.

            const double T_top    = t_eff * (d_j * d_j/(d_j_half * d_j_half)
                - 2.0 * d_j/d_j_half) + t_pole;                         // K at i=im-1
            // ATNEPT_T_INIT_ANCHOR=1: the adiabat through t_1bar_K() at 1 bar, the same at
            // every latitude, isothermal at T_top above the level where it has cooled to T_top.
            const bool anchor = t_init_anchor() != 0;
            const double T_bottom = anchor ? t_bottom_anchor_K() : T_top + gam * L_atm;   // K at i=0

            for(int i = 0; i < im; i++){
                height = get_layer_height(i);                           // km
                t.x[i][j][k] = (T_bottom - gam * height) / t_ref;
                if(anchor) t.x[i][j][k] = std::max(T_bottom - gam * height, T_top) / t_ref;

/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     height[km] = " << height << endl
        << "     t_equator[°C] = " << t_equator - t_ref 
        << "     t_pole[°C] = " << t_pole - t_ref << endl
        << "     gam[K/km] = " << gam
        << "     gam * height[K] = " << gam * height << endl
        << "     t_u[°C] = " << t.x[i][j][k] * t_ref - t_ref
        << "     t_u[K] = " << t.x[i][j][k] * t_ref << endl;
*/

            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_Temperature\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_Temperature ended" << endl;
    return;
}
/*
*
*/
/*
 * ATNEPT_TW_BALANCE=<strength>: THE ZONAL WIND IN THERMAL-WIND BALANCE WITH THE TEMPERATURE.
 * Default 1 with ATNEPT_GEO_BALANCE (the default), else 0. Ported from ATOM_Precipitation (VelocityInitializer.h, balance_thermal_wind),
 * whose note gives the reason it is an initial condition: geostrophic adjustment takes ~1/f,
 * tens of thousands of iterations, and cannot be waited for.
 *
 * theta is the colatitude and w the zonal wind (east positive), so
 *
 *     dw/dz = + g / (a f T) * dT/dtheta,     f = 2 omega cos(theta)
 *
 * The wind at the ANCHOR level is kept (ATNEPT_TW_ANCHOR=<level>; default each column's jet
 * level, where the prescribed |w| is largest; ATOM anchors at the surface) and the shear the model's own
 * temperature implies is integrated away from it, trapezoid in height, at the given strength.
 * Tapered out towards the equator with sin2(lat) / (sin2(lat) + sin2(lat_min)), lat_min =
 * ATNEPT_TW_LATMIN (15 deg), because f -> 0 there; capped at ATNEPT_TW_WMAX (600 m/s). Only w
 * is written (and wn); v and u are left alone.
 */
void cNeptuneModel::balanceThermalWind(){
    static const double s = [](){ const char* e = getenv("ATNEPT_TW_BALANCE"); return e ? atof(e) : (geo_balance() ? 1.0 : 0.0); }();
    if(s == 0.0) return;
    static const double lat_min_deg = [](){ const char* e = getenv("ATNEPT_TW_LATMIN"); return e ? atof(e) : 15.0; }();
    // 600 m/s here, not ATURAN's 400: Neptune's prescribed equatorial jet reaches -437 m/s, and at
    // 400 the cap cut 187 720 cells below their own anchor wind (ATJUP/satchk/giants/NGA).
    static const double w_max = [](){ const char* e = getenv("ATNEPT_TW_WMAX"); return e ? atof(e) : 600.0; }();
    static const int anchor_env = [](){ const char* e = getenv("ATNEPT_TW_ANCHOR"); return e ? atoi(e) : -1; }();
    // the anchor: a fixed level when the knob names one, otherwise each column's JET LEVEL, the
    // level where the prescribed |w| is largest (the initialiser puts the measured wind there)
    const bool jet_anchor = !(anchor_env >= 0 && anchor_env < im);
    tw_anchor_i.assign(jm * km, im - 1);
    const double s2_min = std::pow(std::sin(lat_min_deg * M_PI / 180.0), 2);
    const double dz = L_atm * 1.0e3 * dr;                           // layer thickness [m]
    const double w_lim = w_max / u_0;

    long n_cols = 0, n_capped = 0;
    double worst_dw = 0.0;
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:n_cols,n_capped) reduction(max:worst_dw)
    for(int j = 1; j < jm-1; j++){
        for(int k = 0; k < km; k++){
            const double slat = std::cos(the.z[j]);                 // sin(latitude)
            const double f = 2.0 * omega * slat;
            if(f == 0.0) continue;
            const double wgt = (slat * slat) / (slat * slat + s2_min);
            // shear at a level [(m/s)/m]
            auto shear = [&](int i)->double{
                const double T = t.x[i][j][k] * t_ref;
                if(!(T > 0.0)) return 0.0;
                const double dTdthe = (t.x[i][j+1][k] - t.x[i][j-1][k]) * t_ref / (2.0 * dthe);
                const double a_m = metricRadius(rad.z[i]) * L_atm * 1.0e3;
                const double sh = g * dTdthe / (a_m * f * T);
                return std::isfinite(sh) ? sh : 0.0; };
            int ia = anchor_env;
            if(jet_anchor){
                ia = 0;
                for(int i = 1; i < im; i++) if(std::fabs(w.x[i][j][k]) > std::fabs(w.x[ia][j][k])) ia = i;
            }
            tw_anchor_i[j * km + k] = ia;
            const double w0 = w.x[ia][j][k];
            auto put = [&](int i, double acc){
                double dw = s * wgt * acc;                          // [m/s]
                if(!std::isfinite(dw)) dw = 0.0;
                double wnew = w0 + dw / u_0;
                if(std::fabs(wnew) > w_lim){ wnew = (wnew > 0.0 ? w_lim : -w_lim); n_capped++; }
                if(std::fabs(dw) > worst_dw) worst_dw = std::fabs(dw);
                w.x[i][j][k] = wnew; wn.x[i][j][k] = wnew; };
            double acc = 0.0;
            for(int i = ia + 1; i < im; i++){ acc += 0.5 * (shear(i-1) + shear(i)) * dz; put(i, acc); }
            acc = 0.0;
            for(int i = ia - 1; i >= 0; i--){ acc -= 0.5 * (shear(i+1) + shear(i)) * dz; put(i, acc); }
            n_cols++;
        }
    }
    double wmx = 0.0, wmn = 0.0;
    for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
        wmx = std::max(wmx, w.x[i][j][k]); wmn = std::min(wmn, w.x[i][j][k]); }
    printf("      ATNEPT: thermal-wind balance, strength %.3f, anchor %s (level %d at the equator), lat_min %.1f deg: %ld columns, largest shear added %.2f m/s, w now %.2f .. %.2f m/s, capped cells %ld\n",
           s, jet_anchor ? "each column's jet level" : "fixed level", tw_anchor_i[((jm-1)/2) * km], lat_min_deg, n_cols, worst_dw, wmn * u_0, wmx * u_0, n_capped);
}

void cNeptuneModel::init_PressureDynamic(){
    cout << endl << "      ATNEPT: init_PressureDynamic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ATNEPT_PDYN_INIT (2026-10-12): what the dynamic pressure starts from. 1 is the kinetic energy
    // density of the prescribed wind, which this routine has always written; 0 starts it at zero
    // and leaves it to the pressure solver -- ATSAT's repair of the same day (ATSAT_PDYN_INIT,
    // InitVariables_Sat.cpp has the derivation): p_dyn enters the momentum equations through its
    // gradient, and the radial gradient of that kinetic-energy field is a force nothing balances.
    // ATJUP has always started p_dyn at zero.
    // DEFAULT 0 SINCE 2026-10-12 on the user's word (1 restores, byte-identical to a0eab28).
    // On Neptune the field peaked at 1.53 bar at the equator (r_mix 5.05 kg/m3, w up to 420 m/s).
    // MEASURED, 224 iterations, 4 threads, 1 against 0 (ATJUP/satchk/giants/NP1, NP0):
    //     continuity residuum at iteration 16 / 112 / 224   0.148 / 0.068 / 0.052 -> 0.007 / 0.005 / 0.006
    //     radial wind at 224                                -5.42 .. +2.40 -> -0.11 .. +1.15 m/s
    //     H2O / NH3 / H2S column                            -0.83 / -0.42 / -0.30 % -> +0.15 / +0.06 / +0.05 %
    // v, w and the temperature extrema agree within 3 %, the mean OLR is 0.806 W/m2 in both.
    static const int pdyn_init = [](){ const char* e = getenv("ATNEPT_PDYN_INIT"); return e ? atoi(e) : 2; }();
    // ATNEPT_PDYN_GEO=1: THE BAROTROPIC PRESSURE OF A GEOSTROPHIC START. With the radial forces
    // carried by p_hydro, what is left for p_dyn at the start is the meridional balance: the
    // Coriolis force and the metric term w2 cot(theta)/r on the zonal wind need a meridional
    // pressure gradient, the same at every height because the wind below the jet level is. As
    // rhs_v forms them,
    //     scale_p / rm * dp/dtheta = Coriolis * scale_Cor * 2 omega cos(theta) w + w2 cot(theta) / rm
    // integrated from the north pole with the anchor wind of each column (balanceThermalWind),
    // trapezoid in theta. The centrifugal term of rhs_v is 1e-4 of these and is left out.
    if(pdyn_geo() && metric_split()){
        const double scale_p = 1.0e5 / (r_mix * u_0 * u_0);
        static const double sin_floor = [](){ const char* e = getenv("ATNEPT_SINTHE_FLOOR"); return e ? atof(e) : 0.4; }();
        // ATNEPT_PDYN_GEO_EXACT (2026-10-14; default follows ATNEPT_GEO_BALANCE, 0 = the trapezoid): integrate so that the
        // centred difference rhs_v takes returns the integrand exactly at every integrated row,
        // p[j+1] = p[j-1] + 2 dthe G[j], seeded with the trapezoid's first step. The trapezoid's
        // centred difference is the 1-2-1 mean of G, measured 5e-3 off at the steep jet flank
        // 68-74 N (ATNEPT_VBUDGET, NBR224). Two interleaved chains: a 2-cell mode in p is
        // invisible to the centred difference that reads it.
        static const bool geo_exact = [](){ const char* e = getenv("ATNEPT_PDYN_GEO_EXACT"); return e ? atoi(e) != 0 : geo_balance(); }();
        auto integrate = [&](const std::vector<double>& G, std::vector<double>& P){
            P[0] = 0.0;
            for(int j = 1; j < jm; j++)
                P[j] = (geo_exact && j >= 2) ? P[j-2] + 2.0 * dthe * G[j-1]
                                             : P[j-1] + 0.5 * (G[j-1] + G[j]) * dthe;
        };
        double pmax = 0.0, pmin = 0.0;
        #pragma omp parallel for schedule(static) reduction(max:pmax) reduction(min:pmin)
        for(int k = 0; k < km; k++){
            std::vector<double> G(jm, 0.0), P(jm, 0.0);
            for(int j = 0; j < jm; j++){
                const int ia = ((int)tw_anchor_i.size() == jm * km) ? tw_anchor_i[j * km + k] : im - 1;
                const double rm = metricRadius(rad.z[ia]);
                const double sth = std::max(sin_floor, std::abs(sin(the.z[j])));          // the RHS's floored sine
                const double cth = cos(the.z[j]);
                const double wa = w.x[ia][j][k];
                G[j] = (j == 0 || j == jm-1 || sth < 1.0e-6) ? 0.0
                    : (Coriolis * coriolis_scale() * 2.0 * omega * cth * wa * rm + wa * wa * cth / sth) / scale_p;
            }
            integrate(G, P);
            for(int j = 0; j < jm; j++){
                for(int i = 0; i < im; i++) p_dyn.x[i][j][k] = P[j];
                pmax = std::max(pmax, P[j]); pmin = std::min(pmin, P[j]);
            }
        }
        printf("      ATNEPT: init_PressureDynamic: p_dyn is the barotropic pressure that balances Coriolis and metric force of the zonal wind in the meridional equation; %.6f .. %.6f bar (ATNEPT_PDYN_GEO=0 starts it at zero)\n", pmin, pmax);
        // ATNEPT_PDYN_GEO=2: the rest of the balance, level by level (see pdyn_geo_mode()). rhs_v
        // takes -scale_p dp_dyn/dtheta / rm - dp_hydro/dtheta / rm, so p_baro = Q_i - scale_p p_dyn
        // - p_hydro(start), where Q_i is the same trapezoid of the same integrand as above with the
        // wind and the metric radius of level i. Then the start is balanced at every level as
        // closely as it is at the anchor level (the centred difference of the trapezoid is the
        // 1-2-1 mean of the integrand, the error left at the steep jet flanks).
        if(pdyn_geo_mode() == 2){
            p_baro.assign((size_t)im * jm * km, 0.0);
            computeBuoyancyRefLevel();                          // the buoyancy's base state, if ATNEPT_HYDRO_SPLIT puts it in p_hydro
            computeHydrostaticPressure();                       // p_baro still zero: the metric / Coriolis integral (and the buoyancy with ATNEPT_HYDRO_SPLIT)
            double bmax = 0.0, bmin = 0.0;
            #pragma omp parallel for schedule(static) reduction(max:bmax) reduction(min:bmin)
            for(int k = 0; k < km; k++){
                std::vector<double> G(jm, 0.0), Q(jm, 0.0);
                for(int i = 0; i < im; i++){
                    const double rm = metricRadius(rad.z[i]);
                    for(int j = 0; j < jm; j++){
                        const double sth = std::max(sin_floor, std::abs(sin(the.z[j])));
                        const double cth = cos(the.z[j]);
                        const double wi = w.x[i][j][k];
                        G[j] = (j == 0 || j == jm-1 || sth < 1.0e-6) ? 0.0
                            : Coriolis * coriolis_scale() * 2.0 * omega * cth * wi * rm + wi * wi * cth / sth;
                    }
                    integrate(G, Q);
                    for(int j = 0; j < jm; j++){
                        const double b_ = Q[j] - scale_p * p_dyn.x[i][j][k] - p_hydro.x[i][j][k];
                        p_baro[((size_t)i * jm + j) * km + k] = b_;
                        bmax = std::max(bmax, b_); bmin = std::min(bmin, b_);
                    }
                }
            }
            computeHydrostaticPressure();                       // now with p_baro
            printf("      ATNEPT: init_PressureDynamic: ATNEPT_PDYN_GEO=2, p_baro (the per-level rest, added to p_hydro, frozen) %.6e .. %.6e\n", bmin, bmax);
        }
        return;
    }
    // With ATNEPT_METRIC_SPLIT the metric term is balanced by p_hydro, so a p_dyn that balances it
    // too would be an unopposed force: it starts at zero then.
    if(metric_split() && pdyn_init == 2){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++) p_dyn.x[i][j][k] = 0.0;
        cout << "      ATNEPT: init_PressureDynamic: p_dyn starts at zero, the metric term is carried by p_hydro (ATNEPT_METRIC_SPLIT=0 restores the balanced p_dyn)" << endl;
        return;
    }
    // 2 (2026-10-12): p_dyn starts in balance with the metric term of the radial momentum
    // equation. rhs_u carries + (v2 + w2) / rm, the centrifugal acceleration of the horizontal
    // flow itself, and the prescribed jet makes it large (w = 400 m/s gives 6.5e-3 m/s2 on
    // Neptune). With p_dyn = 0 nothing opposes it: it pushes u upward in every column, and
    // against a rigid lid that is a divergence in the wall cells. The pressure solver would
    // build the balancing gradient itself (the term is in its source), but one sweep per call
    // does not get there. So it is integrated here, d p_dyn / dr = (v2 + w2) / (rm * scale_p),
    // from zero at the bottom, trapezoid in r, with the RHS's own scale_p = 1e5 / (r_mix u_0^2).
    // DEFAULT SINCE 2026-10-12 on the user's word ("repair ... the raised continuity residuum from
    // the rigid lid"). MEASURED with a per-term print of rhs_u (not kept): at iteration 16 the term
    // is the whole of the level-mean tendency of u, everything else 20 times smaller or more; with
    // 2 the pressure gradient cancels it in the interior, a remainder stays in the first and last
    // cell. 224 iterations, rigid lid on, 0 against 2 (ATJUP/satchk/giants/NL224, NB224):
    //     continuity residuum at 16 / 224   0.0046 / 0.0172 -> 0.0011 / 0.0088
    //     radial wind at 224                0.00 .. +0.86 -> -0.05 .. +0.23 m/s
    //     H2O / NH3 / H2S column            +0.075 / +0.031 / +0.028 % -> +0.032 / -0.0001 / +0.0008 %
    // v, w, the temperature extrema, the OLR and the precipitation are unchanged.
    if(pdyn_init == 2){
        const double scale_p = 1.0e5 / (r_mix * u_0 * u_0);
        double pmax = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(max:pmax)
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                double f_prev = 0.0, p_ = 0.0;
                for(int i = 0; i < im; i++){
                    const double f = (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k])
                                   / (metricRadius(rad.z[i]) * scale_p);
                    if(i > 0) p_ += 0.5 * (f_prev + f) * dr;
                    p_dyn.x[i][j][k] = p_; f_prev = f;
                    if(p_ > pmax) pmax = p_;
                }
            }
        }
        printf("      ATNEPT: init_PressureDynamic: p_dyn balances the metric term (v2 + w2)/rm of rhs_u; max %.6f bar (ATNEPT_PDYN_INIT=0 starts it at zero)\n", pmax);
        return;
    }
    if(pdyn_init == 0){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++) p_dyn.x[i][j][k] = 0.0;
        cout << "      ATNEPT: init_PressureDynamic: p_dyn starts at zero (ATNEPT_PDYN_INIT=1 restores the kinetic-energy field)" << endl;
        return;
    }

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                p_dyn.x[i][j][k] = r_mix/2.0 
                    * (pow((u.x[i][j][k] * u_0), 2.0) 
                    + pow((v.x[i][j][k] * u_0), 2.0) 
                    + pow((w.x[i][j][k] * u_0), 2.0))/3.0 * 1e-5; // in bar = 1e5 Pa


/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     height = " << height << endl
        << "     gam = " << gam 
        << "     g = " << g 
        << "     R_ref[J/(g*K)] = " << R_ref 
        << "     p_ref[bar] = " << p_ref 
        << "     exp_pressure = " << exp_pressure << endl
        << "     t_ref[K] = " << t_ref
        << "     t_ref[°C] = " << t_ref - t_ref << endl
        << "     t_u[K] = " << t_u
        << "     t_u[°C] = " << t_u - t_ref
        << "     p_u[bar] = " << p_stat.x[i][j][k] << endl;
*/
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureDynamic\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_PressureDynamic ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::init_PressureStatic(){
    cout << endl << "      ATNEPT: init_PressureStatic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    int j_half = (jm-1)/2;
    double d_j_half = (double)j_half;
    double t_eff = t_pole - t_equator;                                  // non-dimensional
    const double exp_pressure = g/(gam * R_ref);

    // t_equator and t_pole are the TOP (i=im-1) temperatures [K].
    // T_bottom = T_top + gam * L_atm; gam [K/km], L_atm [km] via get_layer_height.
    // P(i,j,k) = p_bottom * (T(i)/T(0))^n,  n = g/(gam*R_ref).
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double d_j = (double)(j);

            const double T_top    = t_eff * (d_j * d_j/(d_j_half * d_j_half)
                - 2.0 * d_j/d_j_half) + t_pole;                         // K at i=im-1
            // ATNEPT_T_INIT_ANCHOR=1: as in init_temperature. r_mix has been scaled to the
            // anchored bottom state by then, so p_bottom below IS p_bottom_bar().
            const bool anchor = t_init_anchor() != 0;
            const double T_bottom = anchor ? t_bottom_anchor_K() : T_top + gam * L_atm;   // K at i=0
            const double p_bottom = r_mix * R_mix * T_bottom * 1e-5;    // bar highest static pressure
            const double h_trop   = (T_bottom - T_top) / gam;           // km, where the adiabat reaches T_top
            const double p_trop   = p_bottom * pow(T_top / T_bottom, exp_pressure);

            for(int i = 0; i < im; i++){
                double height = get_layer_height(i);                    // km
                t.x[i][j][k] = (T_bottom - gam * height) / t_ref;
                if(anchor && height > h_trop){
                    // isothermal layer: hydrostatic with the same gas constant as the exponent
                    t.x[i][j][k] = T_top / t_ref;
                    p_stat.x[i][j][k] = p_trop * exp(-(height - h_trop) * g / (R_ref * T_top));
                    continue;
                }

                // The base is the RATIO T(i)/T(0), not T(i)/t_ref. p_bottom is the pressure at
                // T_bottom, so the exponent must be anchored there too or the whole field is
                // scaled by (T_bottom/t_ref)^n — 607x on Neptune, which put the top of the
                // domain at 15 bar. T_bottom is used rather than t.x[0][j][k] so the expression
                // does not depend on i = 0 having been written earlier in this same loop.
                p_stat.x[i][j][k] = p_bottom
                    * pow(t.x[i][j][k] * t_ref / T_bottom, exp_pressure);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_PressureStatic ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::Forces(){
    cout << endl << "      ATNEPT: Forces" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    double rm, sinthe, costhe, rmsinthe;

    // ATNEPT_BC_PHI_PERIODIC (ATPhys.h): with 1 the seam meridian k = 0 is computed like any other.
    const bool phi_per = ATPhys::phi_periodic<cNeptuneModel>();
    const int  k_first = ATPhys::phi_first(phi_per);

    for(int i = 1; i < im-1; i++){
        rm = metricRadius(rad.z[i]);
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            rmsinthe = rm * sinthe;
                for(int k = k_first; k < km-1; k++){
                    const int kw = ATPhys::phi_wrap(phi_per, k-1, km);
                    const int ke = ATPhys::phi_wrap(phi_per, k+1, km);

// influence of the Coriolis force
                    double Coriolis_rad = - 2.0 * omega * sinthe * w.x[i][j][k];
                    double Coriolis_the = + 2.0 * omega * costhe * w.x[i][j][k];
                    double Coriolis_phi = + 2.0 * omega * (- costhe * v.x[i][j][k] 
                        + sinthe * u.x[i][j][k]);

// influence of the centrifugal force
                    double dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k])/(2.0 * dr);
                    double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k])/(2.0 * dthe);
                    double dpdphi = (p_dyn.x[i][j][ke] - p_dyn.x[i][j][kw])/(2.0 * dphi);

// acting forces
                    CoriolisForce.x[i][j][k] = Coriolis * r_mix 
                        * sqrt((pow(Coriolis_rad, 2) 
                        + pow(Coriolis_the, 2) 
                        + pow(Coriolis_phi, 2))/3.0);

                    CentrifugalForce.x[i][j][k] = centrifugal * r_mix 
                        * omega * omega * rm * (1.0 + fabs(sinthe));

                    BuoyancyForce.x[i][j][k] = buoyancy 
//                        * r_mix * g * (1.0 - (t.x[i][j][k] - 1.0))  //  rho0 * g - rho0 * (t - t0)/t0 * g    for   del_rho << rho0
                        * r_mix * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])  //  rho * g
                        /(r_mix * R_mix * t.x[i][j][k] * t_ref) * 1e5;  // in N/m³

                    PresGradForce.x[i][j][k] = 
                        - sqrt((pow(dpdr, 2) 
                        + pow(dpdthe/rm, 2) 
                        + pow(dpdphi/rmsinthe, 2))/3.0)/L_atm * 1.0e5;

/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     CoriolisForce[kN/m²] = " << CoriolisForce.x[i][j][k] << endl
        << "     CentrifugalForce[N/m³] = " << CentrifugalForce.x[i][j][k] << endl
        << "     BuoyancyForce[N/m³] = " << BuoyancyForce.x[i][j][k] << endl
        << "     PresGradForce[kN/m²] = " << PresGradForce.x[i][j][k] << endl << endl;
*/
            }
        }
    }
/*
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
*/
    cout << "      ATNEPT: Forces ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::init_vapour_cloud_ice(std::string gas, 
    double &c_tropopause, double &coeff_A, double &coeff_B, 
    double &coeff_A_i, double &coeff_B_i, 
    double &t_0, double &t_00, 
    double &ep, double &r, double &m,
    double &C, double &L0, double &R, 
    double &del_alf, double &del_bet,
    Array &c, Array &cloud, Array &ice, Array &cloudiness){

    cout << endl << "      ATNEPT: init_vapour_cloud_ice of "  << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    double t_u = 0.0;
    double p_u = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double magnus = 0.0;
    double r_max_equator = 0.0;
    double r_max_pole = 0.0;
    double r_max_add_equator = 0.0;
    double r_max_add_pole = 0.0;
    double t_add_equator = 0.0;
    double t_add_pole = 0.0;

    if(gas == "CH4"){
        r_max_equator = r_ch4;
        r_max_pole = 0.7 * r_ch4;
        magnus = 1e3;
    }
    if(gas == "H2S"){
        r_max_equator = r_h2s;
        r_max_pole = 0.7 * r_h2s;
        magnus = 1e3;
    }
    if(gas == "H2O"){
        r_max_equator = r_h2o;
        r_max_pole = 0.7 * r_h2o;
        magnus = 1e3;
    }
    if(gas == "NH3"){
        r_max_equator = r_nh3;
        r_max_pole = 0.7 * r_nh3;

        r_max_add_equator = r_nh3_add;
        r_max_add_pole = 0.7 * r_nh3_add;

        t_add_equator = t_0_nh3;
        t_add_pole = 0.7 * t_add_equator;

        magnus = 1e3;
    }

// water vapour distribution decreases approaching tropopause,  
//    at p_stat = 1.0 bar = 1013.0 hPa ... E_Rain = 0.0061 bar = 6,1 hPa (t_u = 0°C) ... E_Rain = 0.00564 bar = 56.4 hPa (t_u = -35°C)

//       initial distribution of  humid atmosphere, cloud and ice

    r_max = std::vector<double>(jm, r_max_pole); // radial location of r_max
    double r_max_eff = r_max_pole - r_max_equator;  // coefficient for the zonal parabolic r_max extention

    r_max_add = std::vector<double>(jm, r_max_add_pole); // radial location of r_max
    double r_max_add_eff = r_max_add_pole - r_max_add_equator;  // coefficient for the zonal parabolic r_max extention

    t_add = std::vector<double>(jm, t_add_pole); // radial location of r_max
    double t_add_eff = t_add_pole - t_add_equator;  // coefficient for the zonal parabolic r_max extention

    double d_j_half = (double)(jm-1)/2.0;

    for(int j = 0; j < jm; j++){
        double d_j = (double)j;
        r_max[j] = r_max_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + r_max_pole;
        r_max_add[j] = r_max_add_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + r_max_add_pole;
        t_add[j] = t_add_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + t_add_pole;
    }

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                t_u = t.x[i][j][k] * t_ref;
                p_u = p_stat.x[i][j][k];

//                E_Rain = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_A, coeff_B);  // saturation species vapour pressure for the liquid phase at t > 0°C in bar
                E_Rain = SaturationAdjustmentNept::saturation_vapour_pressure
                    (t_u, C, L0, R, del_alf, del_bet);

//                q_Rain = ep * E_Rain/(p_u - E_Rain);  // species vapour amount at saturation with species formation in kg/kg
                q_Rain = ep * E_Rain/p_u;  // species vapour amount at saturation with species formation in kg/kg

                c.x[i][j][k] = magnus * r_mix * q_Rain;            // coefficient is arbitrary

                if(c.x[i][j][k] >= r_max[j])  c.x[i][j][k] = r_max[j];

                if((gas == "NH3")&&(p_u >= p_00_nh3))  c.x[i][j][k] = r_max[j];
                if((gas == "H2S")&&(p_u >= p_00_h2s))  c.x[i][j][k] = r_max[j];
                if((gas == "H2O")&&(p_u >= p_0_h2o))   c.x[i][j][k] = r_max[j];
                if((gas == "CH4")&&(p_u >= p_0_ch4))   c.x[i][j][k] = r_max[j];

            }

//                if(t_u <= t_00)  c.x[i][j][k] = 0.0;
/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if((j == 90)&&(k == 180)) cout << endl
                    << gas << "  cloudwater and cloudice    °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k << endl
                    << "   height = " << get_layer_height(i) << endl 
                    << "   cloud_loc = " << cloud_loc[j]
                    << "   cloud_loc_equator = " << cloud_loc_equator 
                    << "   cloud_loc_pole = " << cloud_loc_pole << endl
                    << "   x_cloud = " << x_cloud
                    << "   alfa_s = " << alfa_s << endl
                    << "   Humility_critical[/] = " << SaturationAdjustmentNept::humility_critical(x_cloud, 1.0, 0.8)
                    << "   del_q_ls[kg/kg] = " << del_q_ls << endl
                    << "   t_u[K] = " << t_u
                    << "   t_u[°C] = " << t_u - t_ref
                    << "   p_u[bar] = " << p_u << endl
                    << "   t_0[K] = " << t_0
                    << "   t_00[K] = " << t_00 << endl
                    << "   h_T[/] = " << h_T << endl
                    << "   E_Rain[bar] = " << E_Rain
                    << "   q_Rain[kg/kg] = " << q_Rain << endl
                    << "   c[kg/kg] = " << c.x[i][j][k]
                    << "   cloud[kg/kg] = " << cloud.x[i][j][k]
                    << "   ice[kg/kg] = " << ice.x[i][j][k] << endl
                    << "   cloudiness[/] = " << cloudiness.x[i][j][k] << endl << endl;
*/
        }// end k
    }// end j

/*
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(is_land(SeaMount, i, j, k)){
                    c.x[i][j][k] = 0.0;
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                }
            }
        }
    }
*/
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_vapour_cloud_ice\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_vapour_cloud_ice of "  << gas << endl;

    return;
}
/*
*
*/
