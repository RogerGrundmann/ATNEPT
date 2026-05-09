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

---

## Usage

### Command-line

```bash
./cli/nept neptune/config_atnept.xml
```

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

## Author

Roger Grundmann — roger.grundmann@web.de
