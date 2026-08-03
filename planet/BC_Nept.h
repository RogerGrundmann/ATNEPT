/*
 * Atmosphere General Circulation Modell (ATNEPT)
 * Standalone boundary-condition class for the Neptune model.
 * Declared as friend of cNeptuneModel so it may access all private members
 * through the stored reference.
 *
 * Header-only: all method bodies are inline.
*/

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class cNeptuneModel;
class Array;

using namespace std;

class BC_Nept {
public:

    explicit BC_Nept(cNeptuneModel& model) : m(model) {}

    void bcRadius();
    void bcTheta();
    void bcPhi();
    void initTropopauseLayers();

private:
    cNeptuneModel& m;
};


// -----------------------------------------------------------------------
// Inline implementation
// -----------------------------------------------------------------------
#include "cNeptuneModel.h"
#include "BoundaryConditions.h"
#include "Utils.h"

using namespace AtomUtils;


/*
 * The three boundary passes are the SHARED BoundaryConditions<Planet> now. What used to be here —
 * three field lists walked with a (4/3,-1/3) extrapolation — is the same algorithm the other two
 * models run; the lists moved to cNeptuneModel (BC_Nept.cpp) because they are Neptune's, and the
 * two places this model genuinely differs are named there as bc_margin() = 0 and
 * bc_default_form() = NEUMANN.
 *
 * BC_Nept stays as the name so no call site changes, and because initTropopauseLayers below is
 * ATNEPT's own and has no counterpart in the shared header.
 */
inline void BC_Nept::bcRadius() { BoundaryConditions<cNeptuneModel>(m).bcRadius(); }
inline void BC_Nept::bcTheta()  { BoundaryConditions<cNeptuneModel>(m).bcTheta();  }
inline void BC_Nept::bcPhi()    { BoundaryConditions<cNeptuneModel>(m).bcPhi();    }


inline void BC_Nept::initTropopauseLayers()
{
    const int jm = m.jm;

    m.tropopause_layers = std::vector<double>(jm, m.tropopause_pole);
    cout << endl << "      ATNEPT: init_tropopause_layers" << endl;

    const int i_max  = m.im - 1;
    const int j_max  = jm - 1;
    const int j_half = j_max / 2;
    const double coeff_pole = 285.0;

    for(int j = j_half; j >= 0; j--){
        double x = coeff_pole * (1.0 - (double)(j_half - j) / (double)j_half);
        m.tropopause_layers[j] = AtomUtils::Agnesi(m.tropopause_equator, x);
        m.tropopause_layers[j] = std::round(m.tropopause_layers[j]
            / m.L_atm * (double)i_max);
        m.tropopause_layers[j] = m.tropopause_equator / m.L_atm * (double)i_max;
    }

    for(int j = j_max; j > j_half; j--)
        m.tropopause_layers[j] = m.tropopause_layers[j_max - j];

    cout << "      ATNEPT: init_tropopause_layers ended" << endl;
}
