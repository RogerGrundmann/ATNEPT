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
environment variable to enable. All of them fill diagnostic arrays; none feeds back into the
temperature equation on this model — there is no `ATNEPT_RAD_COUPLING`, unlike ATJUP.

**Radiation** — grey two-stream, shared `Radiation.h`

| Variable | Default | Effect |
|----------|---------|--------|
| `ATNEPT_RADIATION` | 0 | run the solve; fills `Q_rad`, `radiation`, `epsilon` |
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
| `ATNEPT_LOCAL_RHO`, `ATNEPT_METRIC_RADIUS`, `ATNEPT_COSTHE_ABS`, `ATNEPT_PDYN_UNITS` | — | legacy/behaviour switches |

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

1. **There is an unopposed heating excess, and it is the highest-value open item.** Run to 224
   iterations with the radiation diagnostic on, the τ=1 photosphere settles at **205.74 K against a
   T_eff(in) of 59.28 K** — 146 K too warm — emitting **143× the planet's energy budget** and still
   climbing at +0.34 K/iteration. This is the worst row of the four models. Until it is found, **the
   opacity constants cannot be judged against this model at all**.

2. **That number got worse when a real bug was fixed, and the previous one was not better.** Before
   the methane-viscosity correction this model read 2.860× — an artefact of two errors partly
   cancelling. `mue_ch4` held methane's viscosity in *centipoise* as if it were Pa·s, so the
   mass-weighted `mue_mix` came out ~150× too large, and `mue_mix` sets the species diffusivities and
   hence the diffusive-enthalpy sink in `rhs_t`. That sink, ~150× overweighted, was holding the
   column down against the heating excess above.

3. **`ATNEPT_THERMAL_MASSFLUX` is a measurement instrument, not a fix.** Setting it to 0 removes the
   sink entirely and the model runs away harder (176× at 224 iterations), so the term is
   load-bearing even though its form is questionable: it is a flux times a temperature *gradient
   magnitude*, with an absolute value on one component only, so it cannot change sign to oppose a
   runaway.

4. **The static pressure profile was wrong until recently, and old output reflects it.** Any run
   predating the `init_PressureStatic` fix has its entire pressure field inflated by 607×, with the
   top of the domain at 15 bar rather than 0.025 bar. Photosphere numbers from such runs are
   measuring the ceiling of the grid.

5. **The radiation, precipitation and turbulence modules are diagnostic-only here.** They fill their
   arrays and nothing reads them back: there is no `ATNEPT_RAD_COUPLING`, `S_precip_*` reaches no
   RHS, and `ATNEPT_TURB_COUPLING` defaults to 0.

6. **The grey opacity is Jupiter's calibration, not Neptune's.** `C_cia` and `opac_cal` were tuned so
   that Jupiter's photosphere lands at 0.25–0.35 bar. Nothing has recalibrated them here.

---

## Author

Roger Grundmann — roger.grundmann@web.de
