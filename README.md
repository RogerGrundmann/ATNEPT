# ATNEPT

Neptune atmospheric general circulation model based on the numerical framework of the
ATOM climate model. Solves the 3-D Navier-Stokes equations in a spherical shell extending
from the deep H₂O / NH₃ / NH₄SH cloud decks up through the methane condensation layer
that defines Neptune's visible disk.

Despite being the most distant planet in the Solar System and receiving roughly
1/900th of Earth's solar flux, Neptune hosts the fastest sustained winds ever measured
on any planet — equatorial zonal jets reach 580 m/s, retrograde with respect to the
planetary rotation. The driver is internal: Neptune radiates about 2.6 times as much
energy as it absorbs from the Sun, and that interior heat budget is what powers the
storms (Voyager 2's "Great Dark Spot" of 1989, the bright methane "Scooter", and
successive dark spots tracked by Hubble) and the strongest atmospheric convection
known beyond Earth. Modelling this convection-driven circulation in a 3-D framework
is the principal motivation of ATNEPT.

For all relevant data concerning Neptune and its atmosphere the book *Planetary Sciences*
by Imke de Pater and Jack J. Lissauer was indispensable. ATNEPT also serves as the
reference implementation from which CH₄ chemistry was ported to the warmer planet
models (ATJUP, ATSAT).

---

## Physics & Numerics

- **Domain:** spherical shell, 41 × 181 × 361 grid points (r × θ × φ), ~2.7 million cells
- **Vertical extent:** 550 km atmospheric shell (~13.75 km per radial step)
- **Dynamics:** finite-difference discretisation of the 3-D Navier-Stokes equations in spherical coordinates
- **Time integration:** 4th-order Runge-Kutta (inner loop) with a Poisson pressure solver (outer loop)
- **Parallelism:** OpenMP shared-memory threading
- **Thermodynamics:**
  - Temperature initialised as a parabolic pole-to-pole profile; reference temperature 72.5 K
  - Equatorial value 55.5 K, polar value 63.5 K — reflecting Neptune's observed warm-pole
    asymmetry (likely a signature of the same interior-driven convection that powers the jets)
  - Boussinesq buoyancy approximation
  - Clausius-Clapeyron / Sanchez-Lavega SVP formulation for saturation vapour pressures
  - Mixed-phase (liquid + ice) saturation adjustment with iterative convergence
- **Microphysics:** two-category ice scheme adapted from the COSMO weather-forecast model
- **Boundary conditions:** Voyager 2 IRIS temperature/pressure profiles; QuikSCAT and OSCAR
  observational datasets used as zonal-wind templates
- **Planetary constants:** g = 11.1 m/s², Ω = 1.08 × 10⁻⁴ rad/s (≈ 16.11 h sidereal day),
  reference wind speed 250 m/s, mean gas constant R = 3.181 J/(g·K)

---

## Chemical Species

| Species | Phases modelled | Role on Neptune |
|---------|-----------------|-----------------|
| CH₄ | vapour · cloud · ice | Defines the visible cloud deck near the tropopause; methane absorption gives Neptune its blue tint |
| H₂O | vapour · cloud water · cloud ice | Deep tropospheric cloud (~50 bar) |
| NH₃ | vapour · cloud · ice | Mid-troposphere reservoir; sequestered as NH₄SH |
| H₂S | vapour | Required for NH₄SH formation |
| NH₄SH | vapour | Heterogeneous reaction NH₃ + H₂S → NH₄SH at ~230 K |

Pole symmetry of the chemistry fields is checked at iteration 32 by comparing j = 1
against j = jm − 2 in ParaView; theta-derivative scalars are wrapped in `std::abs()`
inside `DiffMassFluxNept` so that pole-antisymmetric gradient components do not break
the symmetry of the resulting flux scalars.

---

## Repository Layout

```
ATNEPT/
├── planet/          # core model (RHS, RK4, thermodynamics, chemistry, I/O)
├── lib/             # array types, config parser, FFT, utilities
├── cli/             # command-line driver (nept)
├── python/          # Cython bindings (pyatnept)
├── neptune/         # run directory (XML config, observational data, output)
│   ├── oscar/       # OSCAR ocean-current dataset (used as zonal-wind template)
│   └── windspeed/   # QuikSCAT surface wind data
├── tinyxml2/        # vendored XML library
├── param.py         # code-generation script (auto-generates parameter files)
└── Makefile
```

---

## Build

**Dependencies:** C++11 compiler with OpenMP support, Python 3, Cython, NumPy.

```bash
# Generate parameter files and build CLI + Python extension
make

# CLI binary only
make nept

# Python extension only
make python

# Clean
make clean
```

The `param.py` script auto-generates several `.inc` / `.pyx` files that parameterise the
model; it runs automatically as part of the build whenever `param.py` itself changes.

### Shared physics headers

Eleven headers in `planet/` are **byte-identical across ATJUP, ATSAT, ATNEPT and ATURAN**:

```
ATPhys.h  BoundaryConditions.h  ConvectiveAdjustment.h  FluxLimiter.h  ParaViewWriter.h
Precipitation.h  PressureSolver.h  Radiation.h  Reporting.h  SaturationAdjustment.h  Turbulence.h
```

There is no submodule and no symlink holding them together — Synology Drive has silently reverted
a working tree once, and a submodule costs friction on every clone. They are plain copies, and
`planet/SHARED.md5` is what makes a divergence loud:

```bash
make check-shared
```

Editing one means: edit it in one repo, copy it to the other three, regenerate its line in **all
four** manifests, and rebuild each.

**`make check-shared` cannot catch everything, and this is the part to read before trusting it.**
It verifies a repo against *its own* manifest, so two repos holding different copies of the same
header both report OK — a state that has already occurred once. The check that does catch it is a
diff between repos, which is why the checksum lines are kept sorted by filename:

```bash
diff ../ATJUP/planet/SHARED.md5 planet/SHARED.md5
```

---

## Optional modules

Every module is **off by default** and a stock run is unaffected by its presence. Set the
environment variable to enable. All of them fill diagnostic arrays. Only one of them can feed back
into the temperature equation — `ATNEPT_RAD_COUPLING`, added last of the four models; precipitation
and turbulence remain diagnostic-only here.

**Radiation** — grey two-stream, shared `Radiation.h`

| Variable | Default | Effect |
|----------|---------|--------|
| `ATNEPT_RADIATION` | 0 | run the solve; fills `Q_rad`, `radiation`, `epsilon` |
| `ATNEPT_RAD_COUPLING` | 0.0 | add `Q_rad` to `rhs_t` as `Q_rad·L_rad/(ρ·cp·u_0·t_ref)`. **1.0 is the physically correct value, not a starting point** — see *Known limitations* for why it is invisible at that setting and what the sweep measured |
| `ATNEPT_SOLAR` | 1 | absorbed shortwave channel (only acts with `ATNEPT_RADIATION`) |
| `ATNEPT_SOLAR_STRENGTH` | 1.0 | scale the absorbed insolation |
| `ATNEPT_SW_TAU_PER_BAR` | 1.0 | move the shortwave absorption level |
| `ATNEPT_CIA_STRENGTH` | 1.0 | scale the H₂/He collision-induced opacity |
| `ATNEPT_OPACITY_STRENGTH` | 1.0 | scale the gas-band and cloud opacity |

**Microphysics and turbulence**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATNEPT_PRECIP` | 0 | precipitation scheme; fills `P_*`, `Q_precip`, `S_precip_*` (diagnostic only here) |
| `ATNEPT_TURB` | 0 | run the closure |
| `ATNEPT_TURB_MODEL` | *param* | override `turb_model` (`k_epsilon`, `k_omega`, `k_omega_SST`) |
| `ATNEPT_TURB_COUPLING` | 0.0 | feed the eddy viscosity into momentum, heat and species diffusion |
| `ATNEPT_CONV_ADJ` | 0 | dry convective adjustment |
| `ATNEPT_SATADJ` | *see code* | mirrored saturation adjustment |
| `ATNEPT_CHEM_ENTHALPY` | *see code* | reaction-enthalpy branch in `thermalmassflux` |

**Numerics and experiment knobs**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATNEPT_THERMAL_MASSFLUX` | 1.0 | scale the diffusive-enthalpy sink in `rhs_t` — see *Known limitations* |
| `ATNEPT_SINTHE_MIN` | 0.0 | env floor on sin θ — **not the value in force**: the integrator uses a hardcoded `sinthe_min = 0.4`, so this accessor is not consulted by default |
| `ATNEPT_PRESS_SOLVER` | *see code* | pressure-solver selection |
| `ATNEPT_STEADY` | 1 | steady-state query in the report |
| `ATNEPT_METRIC_RADIUS` | **24622** | Neptune's mean radius in km, referring the 1/r metric factors to the planet rather than to `rad.z`'s 1..2. **ON by default** — set to `0` for the unshifted metric, which is bit-identical to the pre-flip default |
| `ATNEPT_LOCAL_RHO`, `ATNEPT_COSTHE_ABS`, `ATNEPT_PDYN_UNITS` | — | legacy/behaviour switches |

---

## Diagnostics

Printed at the `checkpoint` cadence:

- **printMinMax** — max/min with location for every prognostic and diagnostic field, including the
  radiation trio (`radiation`, `Q_rad`, `emissivity`), the turbulence fields and the precipitation
  fluxes. All read zero when their module is off, and are printed regardless: an absent array and a
  zero array look the same, and only one of them means the writer works.
- **Equatorial column profile** — `i`, `p[bar]`, `T[K]`, `eps`, `netRad`, `Q_rad` from the model top
  down to the deep boundary. **This model is the reason the table exists.** Its photosphere was
  reported at 17.2287 bar and read as an opacity failure across two commits; the real cause was
  `init_PressureStatic` anchoring `p_bottom` against `T/t_ref` instead of `T/T_bottom`, inflating the
  whole pressure field by 607× and putting the **top of the domain at 15 bar**. One column of
  `p[bar]` shows that at a glance; no column-averaged number can.
- **Photosphere line** (with `ATNEPT_RADIATION=1`) — mean OLR against the input budget, the τ=1 level
  in bar, and the temperature there beside the blackbody flux it implies. Note `T_eff(in)` already
  includes `F_int`, so on a planet with a large internal flux the photosphere should still approach
  it; the internal flux does not excuse an excess.
- **ParaView** — `Radiation`, `Q_rad_mW_m3` and `Emissivity`, plus the turbulence six and the
  precipitation eleven, in all four views (panorama `.vts`; radial, zonal, longitudinal `.vtk`).

---

## Usage

### Command-line — note the TWO arguments

```bash
./cli/nept <path> <config file name>
./cli/nept . config_atnept.xml            # config in the current directory
```

A single combined path fails with `couldn't load config file inside cNeptuneModel`. This is the
reverse of ATOM's `cli/atm`, and the same convention ATSAT and ATURAN use.

```bash
OMP_NUM_THREADS=12 ./cli/nept . config_atnept.xml > run.log 2>&1
```

Runs in this repository's measurements are single-threaded (`OMP_NUM_THREADS=1`) so that results
are bit-reproducible and byte-comparisons between builds mean something.

**This model used to be non-reproducible above one thread, and that is fixed.** The site was
`PressureSolverNept.h`'s Poisson loop — Gauss-Seidel written in place, with
`#pragma omp parallel for collapse(2) schedule(dynamic, 4)` over the very two indices its stencil
reads across, so cell `(i,j,k)` was read by the thread owning `(i+1,j)` or `(i,j+1)` while its
owner was writing it. `schedule(dynamic)` meant the chunk assignment itself varied with timing, so
the same binary at the same thread count differed run to run. It is now serial; see the comment in
that file. Measured at nm=4:

| configuration | 16t run A vs B | 1t vs 16t |
|---|---|---|
| before, `ATNEPT_PRESS_SOLVER=0` (default) | differ, 1 of 7 files | differ, 6 of 7 files |
| after, `ATNEPT_PRESS_SOLVER=0` (default) | **bit-identical** | **bit-identical** |
| `ATNEPT_PRESS_SOLVER=1` (shared red-black) | **bit-identical** | **bit-identical** |

**Nothing was taken back to get this.** One thread ran `collapse(2)` in lexicographic order
already, so the serial loop reproduces the previous 1-thread answer bit-identically, with not one
differing log line — every single-threaded measurement in this file still stands. `computePressure`
was 0.003 s at 16 threads and is 0.03 s serial.

This defect was found in ATURAN and the same check was then run here. ATSAT and ATJUP never had it:
both deleted their per-planet solvers and bind straight to the shared red-black
`PressureSolver<Planet>`, which is byte-identical in all four repositories.

**The log is byte-comparable at any thread count too.** The saturation-adjustment block —
`i_sat`/`j_sat`/`k_sat`, `iter_prec_found`, and the `p_stat`, `T`, `saturation` and per-species
`humid/cloud/ice` values printed with them — was for a while the only thing that still varied.
It is filled under `#pragma omp critical` in `SaturationAdjustmentNept.cpp` and used to record
whichever cell reached the section *last*, so the winner followed thread arrival order. That was
never a race — the section is properly synchronised and writes reporting variables only, and every
output file was bit-identical across it — but it made logs from different thread counts impossible
to diff.

The winner is now chosen by **position** instead: the loop nest is `k`, then `j`, then `i`, so
serial traversal visits `key = (k·jm + j)·im + i` in increasing order and "the last cell found
wins" is exactly "the largest key wins". Taking the maximum reproduces the single-threaded answer
at any thread count. **Which cell is reported did not change** — at 1 thread the log is unchanged
line for line; at 16 threads it now matches it, and two 16-thread runs match each other.

### Python

```python
import sys
sys.path.insert(0, "neptune")
import pyatnept

model = pyatnept.NeptuneModel()
model.load_config("neptune/config_atnept.xml")
model.run()
```

Output is written as VTK / VTS files for visualisation in ParaView (panorama, sphere,
radial, zonal, and longitudinal cross-sections).

---

## Known limitations

None of these stops a run; all of them affect what a result means.

1. **The photosphere is 146 K too warm, and it is mostly redistribution rather than heating.** Run
   to 224 iterations with the radiation diagnostic on, the τ=1 photosphere settles at **205.74 K
   against a T_eff(in) of 59.28 K**, emitting **143× the planet's energy budget** and still climbing.
   This is the worst row of the four models. Until it is found, **the opacity constants cannot be
   judged against this model at all**.

   The column mean says what kind of fault it is:

   | checkpoint | T(i=0) deep | T(i=20) mid | T(i=40) top | column mean |
   |---|---|---|---|---|
   | 1 | 727.75 | 407.19 | 90.09 | 409.34 |
   | 28 | 618.71 | 429.72 | 206.32 | 421.72 |

   The deep loses **109.04 K**, the top gains **116.23 K**, and the mean rises **12.38 K (+3.02 %)**.
   So there *is* a genuine net heat gain — unlike ATURAN, whose mean **falls** 0.9 % over the same
   run and whose fault is purely redistributive — but it is not what dominates: the vertical
   redistribution is roughly **9× larger** than the net gain. Thermal diffusion flattens the initial
   adiabat, nothing anchors the top of the column to the planet's energy budget, and the photosphere
   drifts toward the column mean; a smaller real heat source sits underneath that.

   **Two ice giants with near-identical `rhs_t` giving opposite signs on the column mean is itself
   an open question**, and it is not explained. `ATNEPT_THERMAL_MASSFLUX`'s six-order dominance of
   `rhs_t` (item 3) is the obvious suspect and no more than a suspect. Anyone attacking item 1
   should separate the two components before assuming a single cause: the anchor that is missing
   (`ATNEPT_RAD_COUPLING` now supplies it, and item 8 measures how far it reaches — not far) and the
   +3 % source, which ATURAN does not have at all.

   **The numbers in this item were taken at `3c6d534` and the baseline has moved since.** Measured
   on current HEAD, 224 iterations, radiation on, single-threaded: OLR/in **135.601** and
   T(τ=1) **203.36 K**, against the 143× and 205.74 K quoted above. Of the five commits in between,
   three were verified to leave every output file bit-identical, so the shift belongs to the other
   two — the Coriolis and centrifugal corrections `024c37f` and `e412b1b`. The fault is unchanged in
   kind and the table above has not been re-measured cell by cell; treat its columns as the shape of
   the problem and item 8's `off` row as the current baseline.

2. **That number got worse when a real bug was fixed, and the previous one was not better.** Before
   the methane-viscosity correction this model read 2.860× — an artefact of two errors partly
   cancelling. `mue_ch4` held methane's viscosity in *centipoise* as if it were Pa·s, so the
   mass-weighted `mue_mix` came out ~150× too large, and `mue_mix` sets the species diffusivities and
   hence the diffusive-enthalpy sink in `rhs_t`. That sink, ~150× overweighted, was resisting the
   flattening in item 1. Correcting it unmasked the drift rather than causing it.

3. **`ATNEPT_THERMAL_MASSFLUX` is a measurement instrument, not a fix.** Setting it to 0 removes the
   sink entirely and the model drifts harder (176× at 224 iterations), so the term is load-bearing
   even though its form is questionable: it is a flux times a temperature *gradient magnitude*, with
   an absolute value on one component only, so it cannot change sign to oppose the flattening. It
   dominates `rhs_t` by six orders of magnitude on this model, which is why item 1 names it as the
   suspect for the +3 % net gain that ATURAN does not have.

4. **The static pressure profile was wrong until recently, and old output reflects it.** Any run
   predating the `init_PressureStatic` fix has its entire pressure field inflated by 607×, with the
   top of the domain at 15 bar rather than 0.025 bar. Photosphere numbers from such runs are
   measuring the ceiling of the grid.

5. **The precipitation and turbulence modules are diagnostic-only here.** They fill their arrays and
   nothing reads them back: `S_precip_*` reaches no RHS and `ATNEPT_TURB_COUPLING` defaults to 0.
   Switching either on changes plots, not physics. Radiation is no longer in this list — it can now
   reach `rhs_t` through `ATNEPT_RAD_COUPLING`, which is off by default; see item 8.

6. **The metric radius was corrected and made the default, and results before that commit are not
   comparable with results after it.** `rad.z` runs 1..2, so an unshifted metric put Neptune's
   surface `L_atm` = 550 km from the centre instead of R = 24622 km, making every horizontal
   derivative 45× too large. `metricRadius()` now defaults to the planet's radius;
   `ATNEPT_METRIC_RADIUS=0` restores the old metric bit-identically (verified, 92/92 output files).

   | quantity | `rad.z` metric | corrected |
   |---|---|---|
   | continuity residuum | 1.217334 | **0.050666** |
   | max \|v\| meridional [m/s] | 22.9098 | **1.3930** |
   | max \|u\| radial [m/s] | 80.1531 | **2.6732** |
   | max \|w\| zonal [m/s] | 68.2959 | 68.5485 |
   | T(τ=1) [K] | 205.74 | 203.36 |
   | OLR / input | 143.4 | 135.6 |

   The continuity residuum falls 24× and nothing goes non-finite. The zonal wind is prescribed and
   is untouched, which is the control saying this is the metric and not a blanket damping. The
   meridional wind falls 16.4× rather than the nominal 45× — the response is nonlinear because the
   dynamics rebalance, and proportionality was never the right test. **This is not a fix for item
   1**: OLR/in moves 5 %, because that fault is vertical and this correction is horizontal.

   Two sites still read `rad.z` raw: `Turbulence.h:490` and `:705`. That file is **shared**, so
   changing it costs the four-repo protocol, and `ATNEPT_TURB` defaults off — so they are inert in
   a stock run. Switching the closure on together with the corrected metric means closing that gap
   first.

7. **The grey opacity is Jupiter's calibration, not Neptune's.** `C_cia` and `opac_cal` were tuned so
   that Jupiter's photosphere lands at 0.25–0.35 bar. Nothing has recalibrated them here.

8. **`ATNEPT_RAD_COUPLING` is the anchor item 1 asks for, and on this model it is nowhere near
   enough.** The term is ATJUP's, ported by way of ATURAN's, and it converts the diagnostic flux
   divergence `Q_rad` into a nondimensional temperature tendency,
   `Q_rad·L_rad/(ρ·cp_mix·u_0·t_ref)`, under ATJUP's `rad_t_max = 0.5` limiter. It defaults to 0.

   **The way back is exact.** With the knob unset, 224 iterations, radiation on: **92 of 92 output
   files bit-identical** against `e412b1b`, the commit before the term existed, and the log differs
   only in the `<program name>` line.

   224 iterations, radiation on, single-threaded:

   | coupling | OLR/in | T(τ=1) | τ=1 [bar] | T(i=40) top | T(i=0) deep | raw \|tendency\| max | cells capped |
   |---|---|---|---|---|---|---|---|
   | **0 (off)** | **135.601** | **203.36** | 0.0576 | 199.73 | 609.80 | — | — |
   | 1.0 | 135.600 | 203.36 | 0.0576 | 199.73 | 609.80 | 5.04e−4 | 0 % |
   | 1e3 | 134.060 | 202.76 | 0.0576 | 199.10 | 609.70 | 0.494 | 0 % |
   | 1e4 | 133.036 | 202.35 | 0.0576 | 198.73 | 608.85 | 4.875 | 3.2 % |
   | 3e4 | 138.297 | 204.31 | 0.0576 | 200.65 | 607.59 | 15.233 | 13.3 % |
   | 1e5 | 139.941 | 204.91 | 0.0576 | 201.23 | 606.21 | 51.476 | 86.7 % |

   **At 1.0 — the physically correct value — the term is live but invisible.** It changes all 92
   output files, and moves OLR/in by 0.001 and no temperature in the first two decimals. That was
   predicted before the run from the radiative relaxation time vastly exceeding the step, and
   matches; a *visible* result at 1.0 would have meant a units error. Neptune is colder than Uranus,
   where the same term was invisible at 1.0 too.

   **The last two rows are the limiter, not the term, and the cap was measured rather than
   inferred.** The last two columns come from an instrumented build that records the raw tendency
   before capping, over the same 224 iterations as the rest of the row. The cap bites no cell at all
   at 1.0 and 1e3, 3.2 % of cells at 1e4, 13.3 % at 3e4 and 86.7 % at 1e5, so **the trend reverses
   exactly where capping stops being marginal**: above 1e4 the term redistributes by where the cap
   bites and *warms* the top instead of cooling it.

   **The cap column was first measured at nm=8 and those numbers were wrong.** This item as first
   committed in `c5e7d74` read 1.37e−4 / 0.137 / 1.374 / 4.121 / 13.736 for the raw tendency and
   0 / 0 / 3.5 / 35.2 / 82.1 % for the capped fraction, from 8-iteration runs printed beside a
   224-iteration table. The raw tendency depends on the `Q_rad`/ρ profile, and that profile moves a
   long way over a run — this column's top goes from 90 K to 200 K — so the short runs understated
   it by roughly 3.6× and mis-stated the 3e4 fraction as 35 % against a true 13 %. The physics
   columns were unaffected: an instrumented run reproduces every OLR/in and T(τ=1) in this table
   exactly, which is what says the diagnostic does not perturb the model. **The conclusion is
   unchanged** — 1.0 and 1e3 clean, 1e4 marginal, the reversal where capping becomes substantial —
   but the numbers behind it are these, not those.

   **Read the two clean rows and one nearly-clean one, and they say the term cannot reach the
   answer.** Across 0 → 1e3 → 1e4, four orders of magnitude of coupling, T(τ=1) falls **203.36 →
   202.35 K — about 1 K** — against the ~144 K it would have to fall to meet T_eff(in) = 59.28 K.
   The response also saturates rather than accumulating: the first factor of 1000 buys 0.60 K and
   the next factor of 10 buys 0.41 K. **The term is directionally right and cannot close item 1 at
   this run length**, more decisively than on ATURAN, where the same sweep moved the photosphere a
   few per cent. It is committed as the missing anchor and as a measurement instrument, not as a fix.

   **The cap starts about a decade lower here than on Uranus, and that is now measured on both.**
   The same instrument run against ATURAN at nm=224 gives raw tendencies of 4.39e−5, 0.0438, 0.428,
   1.252 and 4.200 for the same five couplings, capping 0 %, 0 %, 0 %, 1.1 % and 4.7 % of cells. So
   Uranus stays entirely sub-cap up to 1e4 and only reverses at 1e5, where this model is already
   87 % capped at the same setting. Neptune's raw tendency runs **11–12× Uranus's coupling for
   coupling**, which is what a hotter, thinner-topped column does to `Q_rad`/ρ. The two files
   describe different planets, not a disagreement.

---

## Author

Roger Grundmann — roger.grundmann@web.de
