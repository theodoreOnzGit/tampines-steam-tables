# CLAUDE.md

Guidance for Claude Code (and other AI assistants) working in this repository.

## Project overview

TAMPINES Steam Tables is an in-house Rust implementation of the IAPWS-IF97
steam/water property formulation for the **T**hermo-hydraulic **A**rtificial
intelligence **M**ulti-**P**hase **IN**tegrated **E**mulator **S**ystem
(TAMPINES) solver. Unlike the upstream [rust-steam](https://github.com/marciorvneto/rusteam)
library it draws from, this crate uses **dimensioned units** throughout via the
`uom` crate, and incorporates verification tests against the International Steam
Tables (Kretzschmar & Wagner, 2019).

It also provides steam-turbine and converging-diverging nozzle equations,
including choked (critical) two-phase flow, and powers the secondary loop of an
FHR (Fluoride salt-cooled High-temperature Reactor) educational simulator.

License: GPL-3.0 (OpenFOAM-derived algorithms are included; see README).

## Build, test, run

```bash
cargo build                 # build the library
cargo test                  # run all unit/verification tests (~144 test fns)
cargo test <name>           # run a subset by substring match
cargo run --release --example fhr_sim_v2   # FHR educational simulator
```

On Linux, `ndarray-linalg` uses the system OpenBLAS, so you need:

```bash
sudo apt install libopenblas-dev
```

Windows/macOS targets use the static Intel MKL feature instead (see `Cargo.toml`).

## Code layout

Properties are organised by IAPWS-IF97 region under `src/`:

- `region_1_subcooled_liquid/` — region 1 (subcooled liquid)
- `region_2_vapour/` — region 2 (vapour, incl. metastable subregion)
- `region_3_single_phase_plus_supercritical_steam/` — region 3 + supercritical
- `region_4_vap_liq_equilibrium/` — region 4 (saturation line / VLE)
- `region_5_steam_at_800_plus_degc/` — region 5 (ultra-high-temp steam)

Forward equations are `(p,T)` / `(v,T)` flashes. Backward (inverse) equations
live in `backward_eqn_ph_*`, `backward_eqn_ps_*`, `backward_eqn_hs_*`.

Transport and misc properties: `dynamic_viscosity/`, `thermal_conductivity/`,
`surface_tension/`, `dielectric_constant/`.

User-facing entry points are in `interfaces/` — both a functional-programming
API (`(p,T)`, `(p,h)`, `(p,s)`, `(h,s)` flashes) and an object-oriented
`TampinesSteamTableCV` control-volume wrapper. The region-dispatch logic mostly
lives here.

`steam_turbine_equations/` holds nozzle and turbine equations, including the
choked-flow work (see below). `openfoam_algorithms/` contains reference
OpenFOAM solver ports (rhoPimpleFoam, driftFluxFoam, etc.) intended for future
transient two-phase coupling.

## Choked flow (current focus)

`src/steam_turbine_equations/converging_diverging_nozzles/choked_flow/`
implements critical-flow solvers using the Homogeneous Equilibrium Model (HEM):

- `single_phase_basic_choked_flow.rs` — single-phase choked flow.
- `stagnation_point_within_vle_ph_dome_multiphase.rs` — stagnation state inside
  the p-h VLE dome (two-phase).
- `stagnation_point_outside_vle_ph_dome_multiphase.rs` — stagnation state
  outside the dome (subcooled liquid-like, superheated/supercritical).
- `basic_multiphase_equations.rs` — generic multiphase relations (e.g.
  stagnation properties from throat properties).
- `saturation_lookup_table.rs` — precomputed table seeding the bubble/dew-point
  bisection.

Verification tests are under `.../tests/`, validated against:

- Moody (1975), maximum discharge rate of liquid-vapour mixtures — `moody_*`.
- Zaloudek critical mass flux — `zaloudek_*`. NOTE: these reference values are
  graph-read (digitised) HEM curves, not raw experimental data, so keep mass-flux
  (G) tolerances loose.
- Marviken critical flow tests — `marviken_tests.rs`.

### Known sharp edges

- Near the **bubble point**, near-saturated stagnation states must be routed to
  the in-dome solver, not the subcooled one — the dispatcher handles this and it
  is easy to break.
- HEM has documented limitations near the saturation line (see in-code comments
  and `docs/derivation/`); metastable / non-equilibrium effects are not modelled.

## Conventions

- All public property functions take and return `uom` dimensioned quantities —
  do not introduce bare `f64` SI values at API boundaries.
- Match the existing per-region module structure when adding equations
  (`dimensionless_*`, `gamma_*` / `phi_*` derivatives, `intensive_properties.rs`).
- Add a verification test against steam-table or published reference data for any
  new property or flash path; existing tests document expected accuracy bounds.
- The README `# Changelog` is the project's running history — add an entry there
  when bumping the version in `Cargo.toml`.