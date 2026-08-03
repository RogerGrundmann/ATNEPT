// Content moved to BC_Nept.h (header-only standalone class).
// This file is intentionally empty.

/*
 * The four boundary-condition FIELD LISTS, for the shared BoundaryConditions<Planet>.
 *
 * Transcribed verbatim from the arrays inside BC_Nept.h's bcRadius/bcTheta/bcPhi, which already
 * used this shape — a list walked once — so nothing here is a new decision. They are moved out to
 * the model because the shared header asks the planet which fields it has: ATNEPT's lists carry
 * h2s_cloud, h2s_ice and rho_mix, which ATSAT's do not, and that is the whole reason the lists
 * could not live in the shared file.
 *
 * The theta pass is two lists on purpose. Scalars and tracers are extrapolated to the poles;
 * the mass and diffusive fluxes are ZEROED there instead, because they are singular on the axis.
 * ATSAT extrapolates some of what ATNEPT zeroes — a physics disagreement recorded when the
 * boundary conditions were first shared between ATSAT and ATJUP, and preserved here rather than
 * silently resolved.
 */
#include "cNeptuneModel.h"

std::vector<Array*> cNeptuneModel::bc_fields_radius(){
    return {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s, &h2s_cloud, &h2s_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,   &j_nh3,   &j_nh4sh,
        &jT_h2s,  &jT_nh3,  &jT_nh4sh,
        &w_h2s,   &w_nh3,   &w_nh4sh,
        &massflux_h2s,  &massflux_nh3, &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,   &difflux_nh3,  &difflux_nh4sh,
        &thermalmassflux,
        &rho_mix,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
}

std::vector<Array*> cNeptuneModel::bc_fields_theta_extrap(){
    return {
        &t, &u,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s, &h2s_cloud, &h2s_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,   &j_nh3,   &j_nh4sh,
        &jT_h2s,  &jT_nh3,  &jT_nh4sh,
        &w_h2s,   &w_nh3,   &w_nh4sh,
        &thermalmassflux,
        &rho_mix,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
}

std::vector<Array*> cNeptuneModel::bc_fields_theta_zero(){
    return {
        // v and w FIRST, and they are the reason this list is not just the flux fields: BC_Nept
        // zeroed the two tangential velocities at the poles with four hand-written lines outside
        // its zero list, so transcribing only the list left them extrapolated instead. Caught by
        // the acceptance test — 5 of 7 output files differed, and the two that did not (longal at
        // fixed j, and the panorama) are exactly the ones a theta-boundary error cannot reach.
        &v, &w,
        &massflux_h2s, &massflux_nh3, &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,  &difflux_nh3, &difflux_nh4sh
    };
}

std::vector<Array*> cNeptuneModel::bc_fields_phi(){
    return {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s, &h2s_cloud, &h2s_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,   &j_nh3,   &j_nh4sh,
        &jT_h2s,  &jT_nh3,  &jT_nh4sh,
        &w_h2s,   &w_nh3,   &w_nh4sh,
        &massflux_h2s,  &massflux_nh3,  &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,   &difflux_nh3,   &difflux_nh4sh,
        &thermalmassflux,
        &rho_mix,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
}
