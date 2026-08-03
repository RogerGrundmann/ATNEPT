/*
 * Neptune Atmosphere Circulation Model (ATNEPT)
 * Precipitation microphysics — ATNEPT's binding of the SHARED implementation.
 */
#pragma once
#include "Precipitation.h"
class cNeptuneModel;
typedef Precipitation<cNeptuneModel> PrecipitationNept;
