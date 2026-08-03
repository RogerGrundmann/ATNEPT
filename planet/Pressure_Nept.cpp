// Content moved to PressureSolverNept.h (header-only standalone class).
// This file is intentionally empty.

#include "cNeptuneModel.h"

/*
 * Radial-wall conditions for the pressure projection, mirrored from ATSAT's
 * cSaturnModel::prepareProjectionBoundaries. Required by the SHARED PressureSolver.h,
 * which asks the model to state its own wall treatment rather than assuming one.
 *
 * Reached only when ATNEPT_PRESS_SOLVER=1 selects the shared solver; ATNEPT's own
 * PressureSolverNept remains the default and does not call this.
 */
void cNeptuneModel::prepareProjectionBoundaries(bool rigid_lid){
        #pragma omp parallel for
        for(int j = 1; j < jm-1; j++){          // r-direction
            for(int k = 1; k < km-1; k++){
                if(rigid_lid){
                    aux_u.x[0][j][k]      = 0.0;
                    aux_u.x[im-1][j][k] = 0.0;
                } else {
                aux_u.x[0][j][k] = aux_u.x[3][j][k]
                    - 3.0 * aux_u.x[2][j][k] + 3.0 * aux_u.x[1][j][k];
                aux_u.x[im-1][j][k] = aux_u.x[im-4][j][k]
                    - 3.0 * aux_u.x[im-3][j][k] + 3.0 * aux_u.x[im-2][j][k];
                }

                aux_v.x[0][j][k] = aux_v.x[3][j][k]
                    - 3.0 * aux_v.x[2][j][k] + 3.0 * aux_v.x[1][j][k];
                aux_v.x[im-1][j][k] = aux_v.x[im-4][j][k]
                    - 3.0 * aux_v.x[im-3][j][k] + 3.0 * aux_v.x[im-2][j][k];

                aux_w.x[0][j][k] = aux_w.x[3][j][k]
                    - 3.0 * aux_w.x[2][j][k] + 3.0 * aux_w.x[1][j][k];
                aux_w.x[im-1][j][k] = aux_w.x[im-4][j][k]
                    - 3.0 * aux_w.x[im-3][j][k] + 3.0 * aux_w.x[im-2][j][k];

                rhs_u.x[0][j][k] = rhs_u.x[3][j][k]
                    - 3.0 * rhs_u.x[2][j][k] + 3.0 * rhs_u.x[1][j][k];
                rhs_u.x[im-1][j][k] = rhs_u.x[im-4][j][k]
                    - 3.0 * rhs_u.x[im-3][j][k] + 3.0 * rhs_u.x[im-2][j][k];

                rhs_v.x[0][j][k] = rhs_v.x[3][j][k]
                    - 3.0 * rhs_v.x[2][j][k] + 3.0 * rhs_v.x[1][j][k];
                rhs_v.x[im-1][j][k] = rhs_v.x[im-4][j][k]
                    - 3.0 * rhs_v.x[im-3][j][k] + 3.0 * rhs_v.x[im-2][j][k];

                rhs_w.x[0][j][k] = rhs_w.x[3][j][k]
                    - 3.0 * rhs_w.x[2][j][k] + 3.0 * rhs_w.x[1][j][k];
                rhs_w.x[im-1][j][k] = rhs_w.x[im-4][j][k]
                    - 3.0 * rhs_w.x[im-3][j][k] + 3.0 * rhs_w.x[im-2][j][k];
            }
        }

        #pragma omp parallel for
        for(int k = 1; k < km-1; k++){          // theta-direction
            for(int i = 1; i < im-1; i++){
                aux_u.x[i][0][k] = aux_u.x[i][3][k]
                    - 3.0 * aux_u.x[i][2][k] + 3.0 * aux_u.x[i][1][k];
                aux_v.x[i][0][k] = 0.0;
                aux_w.x[i][0][k] = 0.0;

                aux_u.x[i][jm-1][k] = aux_u.x[i][jm-4][k]
                    - 3.0 * aux_u.x[i][jm-3][k] + 3.0 * aux_u.x[i][jm-2][k];
                aux_v.x[i][jm-1][k] = 0.0;
                aux_w.x[i][jm-1][k] = 0.0;

                rhs_u.x[i][0][k] = rhs_u.x[i][3][k]
                    - 3.0 * rhs_u.x[i][2][k] + 3.0 * rhs_u.x[i][1][k];
                rhs_v.x[i][0][k] = 0.0;
                rhs_w.x[i][0][k] = 0.0;

                rhs_u.x[i][jm-1][k] = rhs_u.x[i][jm-4][k]
                    - 3.0 * rhs_u.x[i][jm-3][k] + 3.0 * rhs_u.x[i][jm-2][k];
                rhs_v.x[i][jm-1][k] = 0.0;
                rhs_w.x[i][jm-1][k] = 0.0;
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < im; i++){            // phi-direction
            for(int j = 0; j < jm; j++){
                aux_u.x[i][j][0] = c43 * aux_u.x[i][j][1] - c13 * aux_u.x[i][j][2];
                aux_u.x[i][j][km-1] = c43 * aux_u.x[i][j][km-2] - c13 * aux_u.x[i][j][km-3];
                aux_u.x[i][j][0] = aux_u.x[i][j][km-1] =
                    (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.0;

                aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
                aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
                aux_v.x[i][j][0] = aux_v.x[i][j][km-1] =
                    (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.0;

                aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
                aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
                aux_w.x[i][j][0] = aux_w.x[i][j][km-1] =
                    (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.0;

                rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
                rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
                rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] =
                    (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.0;

                rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
                rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
                rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] =
                    (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.0;

                rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
                rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
                rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] =
                    (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.0;
            }
        }

}

