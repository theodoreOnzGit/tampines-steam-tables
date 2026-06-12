//! Critical-flow solvers for stagnation states that lie OUTSIDE the p-h VLE
//! dome (single phase at the inlet).
//!
//! This module currently covers the **subcooled-liquid / liquid-like** bucket
//! (left of the dome). On isentropic depressurisation such a state stays
//! single-phase liquid down to the bubble point, then flashes into the dome.
//! The superheated-vapour / supercritical buckets (right of the dome) will be
//! added here later.

use uom::ConstZero;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::pressure::megapascal;
use uom::si::pressure::pascal;

use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
use super::basic_multiphase_equations::bubble_point_pressure_from_entropy;

/// Critical pressure & mass flux for a subcooled-liquid / liquid-like
/// stagnation state (OUTSIDE the dome, left side).
///
/// Precondition: `(p0, h0)` is single-phase on the liquid side — subcooled
/// liquid (`p0 < p_c`, `T0 < T_sat(p0)`) or liquid-like (`p0 >= p_c`,
/// `T0 < T_c`). The caller's dispatcher is responsible for routing here.
///
/// Method (energy-balance max-G, HEM — same `G(p)` curve as the in-dome
/// solver, no sound speed involved):
///   along the isentrope `s = s0`,
///     `G(p) = rho(p,s0) * sqrt( 2 * (h0 - h(p,s0)) )`
///
/// * In the single-phase liquid stretch `[p_bubble, p0]`, `rho ~ const` and the
///   enthalpy drop is tiny, so `G` rises monotonically to its liquid-region
///   maximum **at the bubble point**.
/// * Below the bubble point the flow flashes, the HEM sound speed collapses,
///   and `G` develops a (possibly interior) peak in the two-phase region.
///
/// The choke is the global maximum of `G` along the whole isentrope:
///   `G_crit = max( G(p_bubble),  max_{p in [p_min, p_bubble]} G(p) )`
/// * interior two-phase peak wins -> **flashing choke**
/// * bubble-point value wins      -> **bubble-point choke** (strongly subcooled)
///
/// This avoids `mass_flux_ps_eqm_throat` (finite-difference sound speed +
/// bubble-point clamp) entirely; it only needs smooth `h(p,s0)`, `v(p,s0)`.
#[inline]
pub fn get_critical_pressure_and_mass_flux_subcooled_liquid_ph(
    p0: Pressure,
    h0: AvailableEnergy,
) -> (Pressure, MassFlux) {

    let s0 = s_ph_eqm(p0, h0);
    let p_min = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);

    // mass flux from energy conservation at pressure p (along s = s0)
    let g_of_p = |p_pa: f64| -> MassFlux {
        let p = Pressure::new::<pascal>(p_pa);
        let h = h_ps_eqm(p, s0);
        let ke = h0 - h;                       // kinetic energy per unit mass
        if ke < AvailableEnergy::ZERO {
            return MassFlux::ZERO;             // over-expanded guard
        }
        let rho = v_ps_eqm(p, s0).recip();
        rho * (2.0 * ke).sqrt()                // = rho * u
    };

    // bubble point: where the isentrope first reaches saturation (x = 0).
    // The single-phase liquid stretch above it is monotonic in G, so its
    // maximum is exactly this point — one evaluation, no search needed.
    let p_bubble = bubble_point_pressure_from_entropy(s0);
    let g_bubble = g_of_p(p_bubble.get::<pascal>());

    // two-phase maximum: golden-section over [p_min, p_bubble]. G is unimodal
    // here, so this is robust and needs no derivative of the noisy sound speed.
    //
    // Reference:
    //   Price, C. J., & Robertson, B. L. (2012). Golden Section Search.
    //   In Encyclopedia of Engineering Optimization and Heuristics
    //   (pp. 1-4). Singapore: Springer Nature Singapore.
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;     // 0.618...
    let mut a = p_min.get::<pascal>();
    let mut b = p_bubble.get::<pascal>();
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    for _ in 0..100 {
        if (b - a).abs() < 1.0 { break; }      // 1 Pa bracket width
        let gc = g_of_p(c).get::<kilogram_per_square_meter_second>();
        let gd = g_of_p(d).get::<kilogram_per_square_meter_second>();
        if gc > gd {
            b = d;                             // peak is in [a, d]
        } else {
            a = c;                             // peak is in [c, b]
        }
        c = b - gr * (b - a);
        d = a + gr * (b - a);
    }
    let p_two_phase = Pressure::new::<pascal>(0.5 * (a + b));
    let g_two_phase = g_of_p(p_two_phase.get::<pascal>());

    // global maximum along the isentrope = the choke (critical) condition
    if g_two_phase.get::<kilogram_per_square_meter_second>()
        >= g_bubble.get::<kilogram_per_square_meter_second>()
    {
        (p_two_phase, g_two_phase)             // flashing choke
    } else {
        (p_bubble, g_bubble)                   // bubble-point choke
    }
}
