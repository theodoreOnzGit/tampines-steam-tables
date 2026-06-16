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

### Current effort: near-bubble-point HEM artifact

We are trying to solve the **near-bubble-point HEM artifact** that breaks the
Zaloudek VLE critical-pressure / mass-flux tests.

The original combined canary
`zaloudek_*::generic_multiphase_stagnation::quality_0_05_stagnation` is now
`#[ignore]`d. We are splitting the stagnation round-trip by where the stagnation
state lands relative to the VLE dome:

- `subcooled_outside_dome_stagnation.rs` — stagnation outside the dome (left
  side, Region 1 subcooled liquid). Backward-maps each Zaloudek throat to a
  stagnation state, keeps only `ph_flash_region == Region1`, and runs
  `get_critical_pressure_and_mass_flux_subcooled_liquid_ph`. The 20 genuinely-
  subcooled curves (x_t = 0.05 … 1.00) pass.
- in-dome stagnation (two-phase) — counterpart bucket, still being built.

The **active canary** is now
`subcooled_outside_dome_stagnation::quality_bubble_point_subcooled`
(x_t = 1e-4, throats essentially on the saturated-liquid line; intentionally not
`#[ignore]`d while under investigation). It exercises the worst case of the
saturated-liquid-line artifact. The detailed three-failure-mode writeup lives in
the comment block directly above that test; in short, HEM cannot reproduce the
x≈0 choking line in both mass flux and pressure (mass-flux artifact at 5/10 psia,
11–21% choke-pressure error at 15–200 psia, in both solver branches) and a non-
equilibrium / relaxation model would be needed.

The older combined canary swept x_t = 0.05 over a pressure range; first and last
reference points:

- first: p = 5 psia, G = 64.0497 lb/(s·ft²), h0 = 177.3399 Btu/lb
- last:  p = 3000 psia, G = 14016.4977 lb/(s·ft²), h0 = 795.0739 Btu/lb

Diagnosis so far:

- Stagnation reconstruction is fine: `h0_calc ≈ h0_expected` at every point
  (e.g. 343.98 vs 345.81 Btu/lb at 100 psia).
- The forward solver in
  `choked_flow/mod.rs::get_critical_pressure_and_mass_flux_with_stagnation_props`
  locks onto a spurious root near the bubble point. At 100 psia / x≈0.05 it
  converges to p_throat ≈ 860.3 kPa vs the reference 689.5 kPa (+25%), at
  quality ≈ 0.034.
- `g_energy` (energy balance) and `g_hem` (`mass_flux_ps_eqm_throat`,
  finite-difference dv/dP) never truly cross: at the "converged" point
  g_energy ≈ 3092 but g_hem ≈ 5738, so f = g_energy − g_hem ≈ −2646, nowhere
  near zero. The HEM throat mass flux spikes near the saturated-liquid line, so
  the only sign change the bracket finder sees is across that artifact, not a
  physical choke point. Regula falsi then stalls on the discontinuity
  (retained-endpoint problem) and reports the bogus pressure at max_iterations
  instead of failing.
- Pressure-dependent: 5–75 psia stay within the 5% tolerance; 100 psia is the
  first to break — consistent with the `subcooled` test note (11–21%
  choke-pressure error at 15–200 psia). It is the known HEM limitation near the
  saturation line, not a units or reconstruction bug.

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