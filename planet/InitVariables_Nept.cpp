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
            const double T_bottom = T_top + gam * L_atm;                // K at i=0

            for(int i = 0; i < im; i++){
                height = get_layer_height(i);                           // km
                t.x[i][j][k] = (T_bottom - gam * height) / t_ref;

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
void cNeptuneModel::init_PressureDynamic(){
    cout << endl << "      ATNEPT: init_PressureStatic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

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
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_PressureStatic ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::init_PressureStatic(){
    cout << endl << "      ATNEPT: init_PressureDynamic" << endl;

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
            const double T_bottom = T_top + gam * L_atm;                // K at i=0
            const double p_bottom = r_mix * R_mix * T_bottom * 1e-5;    // bar highest static pressure

            for(int i = 0; i < im; i++){
                double height = get_layer_height(i);                    // km
                t.x[i][j][k] = (T_bottom - gam * height) / t_ref;

                p_stat.x[i][j][k] = p_bottom * pow(t.x[i][j][k], exp_pressure);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_PressureDynamic ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::Forces(){
    cout << endl << "      ATNEPT: Forces" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    double rm, sinthe, costhe, rmsinthe;

    for(int i = 1; i < im-1; i++){
        rm = rad.z[i];
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            rmsinthe = rm * sinthe;
                for(int k = 1; k < km-1; k++){

// influence of the Coriolis force
                    double Coriolis_rad = - 2.0 * omega * sinthe * w.x[i][j][k];
                    double Coriolis_the = + 2.0 * omega * costhe * w.x[i][j][k];
                    double Coriolis_phi = + 2.0 * omega * (- costhe * v.x[i][j][k] 
                        + sinthe * u.x[i][j][k]);

// influence of the centrifugal force
                    double dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k])/(2.0 * dr);
                    double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k])/(2.0 * dthe);
                    double dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1])/(2.0 * dphi);

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
