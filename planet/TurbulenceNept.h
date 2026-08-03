/*
 * Neptune Atmosphere Circulation Model (ATNEPT)
 * Turbulence closure — ATNEPT's binding of the SHARED implementation.
 *
 * The closure lives in Turbulence.h, byte-identical to ATSAT's and ATJUP's copies.
 */

#pragma once

#include "Turbulence.h"

class cNeptuneModel;

typedef Turbulence<cNeptuneModel> TurbulenceNept;
