#include "cNeptuneModel.h"
#include "Reporting.h"

using namespace std;

namespace{
    string heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
    string heading_2 = " results based on three dimensional considerations of the problem";
    string level = "km";

    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    };
    HemisphereCoords convert_coords(double lon, double lat){
        HemisphereCoords ret;
        if(lat > 90){
            ret.lat = lat - 90;
            ret.north_or_south = "°S";
        }else{
            ret.lat = 90 - lat;
            ret.north_or_south = "°N";
        }
        if(lon > 180){
            ret.lon = 360 - lon;
            ret.east_or_west = "°W";
        }else{
            ret.lon = lon;
            ret.east_or_west = "°E";
        }
        return ret;
    }
}
/*
*
*/
void cNeptuneModel::printMinMax(){

    cout << endl << endl << " Courant time step   dt = " << dt << endl << endl;

    cout << endl << endl << " Temperatures " << endl;
    searchMinMax_3D(" max 3D temperature ", " min 3D temperature ", 
//        "  °C", t, 273.15, [](double i)->double{return i - 273.15;}, true);
        "  °C", t, t_ref, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max 3D thermalflux ", " min 3D thermalflux ", "m/s", thermalmassflux, 1.0);
    cout << endl;

    cout << endl << " Velocities " << endl;
    searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", "m/s", w, u_0);
    cout << endl;

    cout << endl << " Pressures and Mixture-Density" << endl;
    searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", "mbar", p_dyn, 1e3);
    searchMinMax_3D(" max 3D pressure static ", " min 3D pressure static ", "bar", p_stat, 1.0);
    searchMinMax_3D(" max 3D rho_mix  ", " min 3D rho_mix  ", "kg/m³", rho_mix, 1.0);
    cout << endl;

    cout << endl << " Water " << endl;
    searchMinMax_3D(" max 3D h2o ",  " min 3D h2o ", " kg/m3", h2o, 1.0);
    searchMinMax_3D(" max 3D h2o_cloud ", " min 3D h2o_cloud ", " kg/m3", h2o_cloud, 1.0);
    searchMinMax_3D(" max 3D h2o_ice ", " min 3D h2o_ice ", " kg/m3", h2o_ice, 1.0);
    cout << endl;

    cout << endl << " Methane " << endl;
    searchMinMax_3D(" max 3D ch4 ",  " min 3D ch4 ", " kg/m3", ch4, 1.0);
    searchMinMax_3D(" max 3D ch4_cloud ", " min 3D ch4_cloud ", " kg/m3", ch4_cloud, 1.0);
    searchMinMax_3D(" max 3D ch4_ice ", " min 3D ch4_ice ", " kg/m3", ch4_ice, 1.0);
    cout << endl;

    cout << endl << " Hydrogen Sulfide " << endl;
    searchMinMax_3D(" max 3D h2s ",  " min 3D h2s ", " kg/m3", h2s, 1.0);
    searchMinMax_3D(" max 3D h2s_cloud ", " min 3D h2s_cloud ", " kg/m3", h2s_cloud, 1.0);
    searchMinMax_3D(" max 3D h2s_ice ", " min 3D h2s_ice ", " kg/m3", h2s_ice, 1.0);
    searchMinMax_3D(" max 3D w_h2s ", " min 3D w_h2s ", " kg/m3", w_h2s, 1.0);
    searchMinMax_3D(" max 3D j_h2s ", " min 3D j_h2s ", " kg/m3", j_h2s, 1.0);
    searchMinMax_3D(" max 3D jT_h2s ", " min 3D jT_h2s ", " kg/m3", jT_h2s, 1.0);
    searchMinMax_3D(" max 3D massflux_h2s ", " min 3D massflux_h2s ", " kg/m3", massflux_h2s, 1.0);
    searchMinMax_3D(" max 3D diff_h2s ", " min 3D diff_h2s ", " kg/m3", difflux_h2s, 1.0);
    cout << endl;

    cout << endl << " Ammonia " << endl;
    searchMinMax_3D(" max 3D nh3 ",  " min 3D nh3 ", " kg/m3", nh3, 1.0);
    searchMinMax_3D(" max 3D nh3_cloud ", " min 3D nh3_cloud ", " kg/m3", nh3_cloud, 1.0);
    searchMinMax_3D(" max 3D nh3_ice ", " min 3D nh3_ice ", " kg/m3", nh3_ice, 1.0);
    searchMinMax_3D(" max 3D w_nh3 ", " min 3D w_nh3 ", " kg/m3", w_nh3, 1.0);
    searchMinMax_3D(" max 3D j_nh3 ", " min 3D j_nh3 ", " kg/m3", j_nh3, 1.0);
    searchMinMax_3D(" max 3D jT_nh3 ", " min 3D jT_nh3 ", " kg/m3", jT_nh3, 1.0);
    searchMinMax_3D(" max 3D massflux_nh3 ", " min 3D massflux_nh3 ", " kg/m3", massflux_nh3, 1.0);
    searchMinMax_3D(" max 3D diff_nh3 ", " min 3D diff_nh3 ", " kg/m3", difflux_nh3, 1.0);
    cout << endl;

    cout << endl << " Ammonia Hydrosufide " << endl;
    searchMinMax_3D(" max 3D nh4sh ",  " min 3D nh4sh ", " mg/m3", nh4sh, 1e6);
    searchMinMax_3D(" max 3D w_nh4sh ", " min 3D w_nh4sh ", " mg/m3", w_nh4sh, 1e6);
    searchMinMax_3D(" max 3D j_nh4sh ", " min 3D j_nh4sh ", " mg/m3", j_nh4sh, 1e6);
    searchMinMax_3D(" max 3D jT_nh4sh ", " min 3D jT_nh4sh ", " mg/m3", jT_nh4sh, 1e6);
    searchMinMax_3D(" max 3D massflux_nh4sh ", " min 3D massflux_nh4sh ", " mg/m3", massflux_nh4sh, 1e6);
    searchMinMax_3D(" max 3D diff_nh4sh ", " min 3D diff_nh4sh ", " mg/m3", difflux_nh4sh, 1e6);
    cout << endl;

    cout << endl << " Forces " << endl;
    searchMinMax_3D(" max 3D Coriolis force ", " min 3D Coriolis force ", "mN/m³", CoriolisForce, 1e3);
    searchMinMax_3D(" max 3D centrifugal force ", " min 3D centrifugal force ", " mN/m³", CentrifugalForce, 1e3);
    searchMinMax_3D(" max 3D buoyancy force ", " min 3D buoyancy force ", "kN/m³", BuoyancyForce, 1e-3);
    searchMinMax_3D(" max 3D presgrad force ", " min 3D presgrad force ", "N/m³", PresGradForce, 1.0);
    cout << endl;

    cout << endl << " Energies " << endl;
    searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", "W/m³", Q_Sensible, 1.0);
    searchMinMax_3D(" max 3D latent heat ", " min 3D latent heat ", "W/m³", Q_Latent, 1.0);
    cout << endl << endl;
}
/*
*
*/
/*
*
*/
// Forwarders to the SHARED Reporting<cNeptuneModel>. The bodies used to be here in full; the
// machinery is identical in all three models. See Reporting.h.
void cNeptuneModel::searchMinMax_3D(string name_maxValue, string name_minValue,
    string name_unitValue, Array &value_D, double coeff,
    std::function< double(double) > lambda, bool print_heading){
    Reporting<cNeptuneModel>(*this).searchMinMax_3D(name_maxValue, name_minValue,
        name_unitValue, value_D, coeff, lambda, print_heading);
}
/*
*
*/
void cNeptuneModel::searchMinMax_2D(string name_maxValue, string name_minValue,
    string name_unitValue, Array_2D &value, double coeff){
    Reporting<cNeptuneModel>(*this).searchMinMax_2D(name_maxValue, name_minValue,
        name_unitValue, value, coeff);
}

double cNeptuneModel::out_maxValue() const{
    return maxValue;
}
/*
*
*/
double cNeptuneModel::out_minValue() const{
    return minValue;
}
/*
*
*/
void cNeptuneModel::print_welcome_msg(){
    if(verbose){
        cout << endl << endl << endl;
        cout << "***** Atmosphere Neptune General Circulation Model(ATNEPT) applied to laminar flow" << endl;
        cout << "***** program for the computation of Neptune-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 6 additional transport equations to describe the water vapour, cloud water, cloud ice and nh3 vapour, nh3 cloud and nh3 ice" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonally constant" << endl;
        cout << "***** water and nh3 vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code(German Weather Forecast)" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)" << endl << endl;
        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_welcome_msg_printed = true;
    }
    return;
}
/*
*
*/
void cNeptuneModel::initMsg(){
    cout << "  present state of the computation " << endl << "  current number of iterations nm = " << nm << "    n = " << n << endl;
    return;
}
/*
*
*/
void cNeptuneModel::print_final_msg(){
    cout << endl << "***** end of the NeptuneAtmosphere General Circulation Modell(ATNEPT) *****" << endl << endl;
    if(n == nm)   cout <<  "***** number of artificial time steps      n = " << n << ", end of program reached because of limit of maximum artificial time steps ***** \n\n" << endl;
}
/*
*
*/
/*
*
*/
/*
 * The four defects ATNEPT's own copy carried, all of which the shared version does not:
 *   - the continuity residual tested fabs(residuum) but stored the SIGNED value, so one negative
 *     residual made every later cell compare true and the reported location became the last cell
 *     scanned;
 *   - min_nh4sh/max_nh4sh were never initialised, so case 13 read the stack;
 *   - the pressure case was commented out with its `break` left OUTSIDE the comment, so it fell
 *     through and printed case 1's continuity residual under the label "dp: pressure Poisson
 *     equation" — p_dynn did not exist to difference against, which is why it was commented out
 *     in the first place. p_dynn is a real array now, maintained by restoreVar with the other
 *     n-copies;
 *   - i_loc_level multiplied by 1.e-3, treating L_atm as metres, so every level printed 0km.
 * None of this had ever been seen, because nothing called the routine.
 */
void cNeptuneModel::steadyQuery(){
    Reporting<cNeptuneModel>(*this).steadyQuery();
}
