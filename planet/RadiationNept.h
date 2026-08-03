/*
 * Neptune Atmosphere Circulation Model (ATNEPT)
 * Grey multi-layer radiation — ATNEPT's binding of the SHARED implementation.
 *
 * The scheme lives in Radiation.h, byte-identical to ATSAT's and ATJUP's copies. Neptune's own
 * numbers are in cNeptuneModel.h, where measured properties of the planet belong.
 */

#pragma once

#include "Radiation.h"

class cNeptuneModel;

typedef Radiation<cNeptuneModel> RadiationNept;
