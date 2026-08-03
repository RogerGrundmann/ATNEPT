/*
 * Neptune Atmosphere Circulation Model (ATNEPT)
 * Dry convective adjustment — ATNEPT's binding of the SHARED implementation.
 *
 * The algorithm, the reasoning and the knobs live in ConvectiveAdjustment.h, which is
 * byte-identical to ATSAT's and ATJUP's copies and knows nothing about which planet it runs on.
 * This file exists only so the call sites keep their familiar name.
 */

#pragma once

#include "ConvectiveAdjustment.h"

class cNeptuneModel;

typedef ConvectiveAdjustment<cNeptuneModel> ConvectiveAdjustmentNept;
