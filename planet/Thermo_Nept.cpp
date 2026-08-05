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

void cNeptuneModel::TropopauseLocation(){
//    cout << endl << "      ATNEPT: TropopauseLocation" << endl;

// parabolic tropopause location distribution from pole to pole assumed
    im_tropopause = std::vector<int>(jm, 0);
    int j_half = (jm-1)/2;
    double d_j_half = (double)j_half;
    double trop_u2_eff = (double)(i_beg - i_max);
//    double trop_u2_eff = (double)(i_beg_trop - i_max_trop);
// computation of the tropopause from pole to pole
    double d_j = 0.0;

//    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        d_j = (double)j;
        im_tropopause[j] = (int)((trop_u2_eff * (d_j * d_j/(d_j_half * d_j_half) 
            - 2.0 * d_j/d_j_half)) + (double)i_beg);
//            - 2.0 * d_j/d_j_half)) + (double)i_beg_trop);
// cout << "   j = " << j << "   im_tropopause[j] = " << im_tropopause[j] << endl;
    }

//    cout << "      ATNEPT: TropopauseLocation ended" << endl;
    return;
}
/*
*
*/
void cNeptuneModel::Latent_Heat(){
    cout << endl << "      ATNEPT: Latent_Heat" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    double Latency_Ice = 0.0; 

    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Q_Latent.x[0][j][k] = 0.0;
            Q_Sensible.x[0][j][k] = 0.0;  // sensible heat in [W/m2] from energy transport equation
        }
    }

    for(int j = 1; j < jm-1; j++){
        double sinthe = sin(the.z[j]);
        for(int k = 1; k < km-1; k++){
            for(int i = im-2; i >= 1; i--){
                double rm = metricRadius(rad.z[i]);
                double rmsinthe = rm * sinthe;

                double t_u = t.x[i][j][k] * t_ref;
                double p_u = p_stat.x[i][j][k];

                double E_Rain = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_h2o_A, coeff_h2o_B);  // saturation h2o vapour pressure for the water phase at t > 0°C in hPa
                double E_Ice = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_h2o_A_i, coeff_h2o_B_i);  // saturation h2o vapour pressure for the ice phase in hPa
                double q_Rain = ep_h2o * E_Rain /(p_u - E_Rain);  // h2o vapour amount at saturation with water formation in kg/kg
                double q_Ice = ep_h2o * E_Ice /(p_u - E_Ice);  // h2o vapour amount at saturation with ice formation in kg/kg

                double E_Rain_h2s = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_h2s_A, coeff_h2s_B);  // saturation h2s vapour pressure for the water phase at t > 0°C in hPa
                double E_Ice_h2s = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_h2s_A_i, coeff_h2s_B_i);  // saturation h2s vapour pressure for the ice phase in hPa
                double q_Rain_h2s = ep_h2s * E_Rain_h2s/(p_u - E_Rain_h2s);  // h2s vapour amount at saturation with water formation in kg/kg
                double q_Ice_h2s = ep_h2s * E_Ice_h2s/(p_u - E_Ice_h2s);  // h2s vapour amount at saturation with ice formation in kg/kg

                double E_Rain_nh3 = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_nh3_A, coeff_nh3_B);  // saturation nh3 vapour pressure for the water phase at t > 0°C in hPa
                double E_Ice_nh3 = 1e3 * SaturationAdjustmentNept::clausius_clapeyron(t_u, coeff_nh3_A_i, coeff_nh3_B_i);  // saturation nh3 vapour pressure for the ice phase in hPa
                double q_Rain_nh3 = ep_nh3 * E_Rain_nh3/(p_u - E_Rain_nh3);  // nh3 vapour amount at saturation with water formation in kg/kg
                double q_Ice_nh3 = ep_nh3 * E_Ice_nh3/(p_u - E_Ice_nh3);  // nh3 vapour amount at saturation with ice formation in kg/kg

                double u_av = 0.5 * (u.x[i+1][j][k] + u.x[i-1][j][k]);
                double v_av = 0.5 * (v.x[i+1][j][k] + v.x[i-1][j][k]);
                double w_av = 0.5 * (w.x[i+1][j][k] + w.x[i-1][j][k]);

                double velocity_av =
                        - sqrt((pow(u_av, 2) 
                        + pow(v_av, 2) 
                        + pow(w_av, 2))/3.0);

                double dtdr = (t.x[i+1][j][k] - t.x[i-1][j][k])/(2.0 * dr);
                double dtdthe = (t.x[i][j+1][k] - t.x[i][j-1][k])/(2.0 * rm * dthe);
                double dtdphi = (t.x[i][j][k+1] - t.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dh2odr = (h2o.x[i+1][j][k] - h2o.x[i-1][j][k])/(2.0 * dr);
                double dh2odthe = (h2o.x[i][j+1][k] - h2o.x[i][j-1][k])/(2.0 * rm * dthe);
                double dh2odphi = (h2o.x[i][j][k+1] - h2o.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dh2sdr = (h2s.x[i+1][j][k] - h2s.x[i-1][j][k])/(2.0 * dr);
                double dh2sdthe = (h2s.x[i][j+1][k] - h2s.x[i][j-1][k])/(2.0 * rm * dthe);
                double dh2sdphi = (h2s.x[i][j][k+1] - h2s.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dnh3dr = (nh3.x[i+1][j][k] - nh3.x[i-1][j][k])/(2.0 * dr);
                double dnh3dthe = (nh3.x[i][j+1][k] - nh3.x[i][j-1][k])/(2.0 * rm * dthe);
                double dnh3dphi = (nh3.x[i][j][k+1] - nh3.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dtemp = dtdr + dtdthe + dtdphi;
                double dh2o = dh2odr + dh2odthe + dh2odphi;
                double dh2s = dh2sdr + dh2sdthe + dh2sdphi;
                double dnh3 = dnh3dr + dnh3dthe + dnh3dphi;


                if(h2o.x[i][j][k] >= q_Rain)  
                    Q_Latent.x[i][j][k] = lv_h2o * velocity_av * dh2o/(L_atm * L_atm) * qheat_fix();
                else  Q_Latent.x[i][j][k] = 0.0;

                if(h2o.x[i][j][k] >= q_Ice)  
                    Latency_Ice = ls_h2o * velocity_av * dnh3/(L_atm * L_atm) * qheat_fix();
                else  Latency_Ice = 0.0;



                if(h2s.x[i][j][k] >= q_Rain_h2s)  
                    Q_Latent.x[i][j][k] = Q_Latent.x[i][j][k] + lv_h2s 
                        * velocity_av * dh2s/(L_atm * L_atm) * qheat_fix();
                else  Q_Latent.x[i][j][k] = 0.0;

                if(h2s.x[i][j][k] >= q_Ice_h2s)  
                    Latency_Ice = Latency_Ice + ls_h2s * velocity_av * dh2s/(L_atm * L_atm) * qheat_fix();
                else  Latency_Ice = 0.0;



                if(nh3.x[i][j][k] >= q_Rain_nh3)  
                    Q_Latent.x[i][j][k] = Q_Latent.x[i][j][k] + lv_nh3 
                        * velocity_av * dnh3/(L_atm * L_atm) * qheat_fix();
                else  Q_Latent.x[i][j][k] = 0.0;

                if(nh3.x[i][j][k] >= q_Ice_nh3)  
                    Latency_Ice = Latency_Ice + ls_nh3 * velocity_av * dnh3/(L_atm * L_atm) * qheat_fix();
                else  Latency_Ice = 0.0;



                Q_Latent.x[i][j][k] = Q_Latent.x[i][j][k] + Latency_Ice;  // latent heat in [W/m³] from energy transport equation



                Q_Sensible.x[i][j][k] = r_mix * cp_mix 
                    * velocity_av * dtemp * t_ref/(L_atm * L_atm) * qheat_fix();  // sensible heat in [W/m³] from energy transport equation


                if(SeaMount.x[i][j][k] == 1.0){
                    Q_Latent.x[i][j][k] = 0.0;
                    Q_Sensible.x[i][j][k] = 0.0;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Latent_Heat\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: Latent_Heat ended" << endl;
    return;
    }
/*
*
*/
