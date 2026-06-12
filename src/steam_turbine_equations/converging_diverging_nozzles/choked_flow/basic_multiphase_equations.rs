use uom::si::f64::*;
use uom::si::pressure::megapascal;
use uom::si::ratio::ratio;
use super::saturation_lookup_table::bubble_point_bracket;
use super::saturation_lookup_table::dew_point_bracket;
use crate::constants::p_crit_water;
use crate::constants::s_crit_water;
use crate::interfaces::functional_programming::hs_flash_eqm::p_hs_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::interfaces::functional_programming::pt_flash_eqm::s_tp_eqm_two_phase;
use crate::prelude::functional_programming::ps_flash_eqm::mass_flux_ps_eqm_throat;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
/// Given throat conditions (p_t, s_t), compute the critical mass flux
/// and back-calculate the stagnation conditions (p_0, h_0)
///
/// This is the inverse of the usual approach — instead of finding
/// the throat from stagnation conditions, we fix the throat and
/// recover the stagnation state.
///
/// From energy conservation (isentropic):
/// h_0 = h_t + G*² / (2 * rho_t²)
///      = h_t + v_t² * G*² / 2
///
/// Entropy is conserved: s_0 = s_t
/// Stagnation pressure recovered via p_hs_eqm(h_0, s_0)
///
/// Reference: Saha (1978) NUREG/CR-0417, eq. 10
///            Moody (1975) NEDO-21052
///
/// Note that this uses the homogeneous equilibrium model.
/// This was validated using Zaloudek's data
#[inline]
pub fn get_stagnation_conditions_from_throat_ps(
    p_t: Pressure,
    s_t: SpecificHeatCapacity,
) -> (Pressure, AvailableEnergy, MassFlux) {

    // critical mass flux at throat
    let g_crit = mass_flux_ps_eqm_throat(p_t, s_t);

    // throat specific volume and enthalpy
    let v_t = v_ps_eqm(p_t, s_t);
    let h_t = h_ps_eqm(p_t, s_t);

    // stagnation enthalpy from energy conservation:
    // h_0 = h_t + 0.5 * u_t²
    // u_t = G* * v_t  (u = G/rho = G*v)
    let u_t: Velocity = g_crit * v_t;
    let h_0 = h_t + 0.5 * u_t * u_t;

    // entropy conserved along isentrope
    let s_0 = s_t;

    // recover stagnation pressure from (h_0, s_0)
    let p_0 = p_hs_eqm(h_0, s_0);

    (p_0, h_0, g_crit)
}

/// Same as above but takes throat (p_t, h_t) as input
/// converts h_t to s_t internally
/// This was validated using Zaloudek's data
#[inline]
pub fn get_stagnation_conditions_from_throat_ph(
    p_t: Pressure,
    h_t: AvailableEnergy,
) -> (Pressure, AvailableEnergy, MassFlux) {
    let s_t = s_ph_eqm(p_t, h_t);
    get_stagnation_conditions_from_throat_ps(p_t, s_t)
}

/// Bubble-point pressure along an isentrope `s = s0`.
///
/// Returns the pressure `p_bubble` at which the saturated-liquid entropy
/// equals `s0` — i.e. the pressure where an isentropic depressurisation of a
/// subcooled / liquid-like state first reaches saturation (x = 0, flashing
/// inception).
///
/// The saturated-liquid entropy
///   `s_f(p) = s_tp_eqm_two_phase(T_sat(p), p, 0.0)`
/// is monotonically increasing in `p` (from ~0 at the triple point up to
/// `s_crit` at the critical point), so the root `s_f(p_bubble) = s0` is unique
/// and recovered by bisection. This automatically handles the Region-3 cap
/// (16.529-22.064 MPa), where the saturated-liquid properties come from the
/// Region 3 EOS.
///
/// Precondition: `s0` lies on the liquid side of the dome, i.e.
/// `s_f(p_triple) <= s0 <= s_crit`. At (or above) the critical entropy the
/// bubble line meets the critical point, so `p_crit` is returned directly.
/// Below the triple-point saturated-liquid entropy it clamps to `p_min`.
#[inline]
pub fn bubble_point_pressure_from_entropy(s0: SpecificHeatCapacity) -> Pressure {
    let p_min = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let p_crit = p_crit_water();

    // at / above the critical entropy the bubble line terminates at the
    // critical point — return the critical pressure directly
    if s0 >= s_crit_water() { return p_crit; }

    // saturated-liquid entropy at pressure p
    let s_f = |p: Pressure| -> SpecificHeatCapacity {
        s_tp_eqm_two_phase(sat_temp_4(p), p, 0.0)
    };

    // below the triple-point saturated-liquid entropy, clamp to p_min
    if s0 <= s_f(p_min) { return p_min; }

    // seed a tight bracket from the saturation lookup table, then refine by
    // bisection within it (s_f is monotonically increasing in p)
    let (mut p_lo, mut p_hi) = bubble_point_bracket(s0);
    for _ in 0..40 {
        let p_mid = 0.5 * (p_lo + p_hi);
        if ((p_hi - p_lo) / p_mid).get::<ratio>() < 1e-9 { break; }
        if s_f(p_mid) < s0 {
            p_lo = p_mid;
        } else {
            p_hi = p_mid;
        }
    }
    0.5 * (p_lo + p_hi)
}

/// Same as [`bubble_point_pressure_from_entropy`] but takes a `(p, h)`
/// stagnation state and uses its entropy. Convenient for reading subcooled /
/// liquid-like points straight off a p-h diagram.
#[inline]
pub fn bubble_point_pressure_ph(p0: Pressure, h0: AvailableEnergy) -> Pressure {
    let s0 = s_ph_eqm(p0, h0);
    bubble_point_pressure_from_entropy(s0)
}

/// Dew-point pressure along an isentrope `s = s0`.
///
/// Returns the pressure `p_dew` at which the saturated-vapour entropy equals
/// `s0` — i.e. the pressure where an isentropic depressurisation of a
/// superheated-vapour / supercritical state first reaches saturation (x = 1,
/// condensation inception). This is the vapour-side analogue of
/// [`bubble_point_pressure_from_entropy`].
///
/// The saturated-vapour entropy
///   `s_g(p) = s_tp_eqm_two_phase(T_sat(p), p, 1.0)`
/// is monotonically *decreasing* in `p` (from large values near the triple
/// point down to `s_crit` at the critical point), so the root
/// `s_g(p_dew) = s0` is unique and recovered by bisection. This handles the
/// Region-3 cap automatically.
///
/// Precondition: `s0` lies on the vapour side of the dome, i.e.
/// `s_crit <= s0 <= s_g(p_triple)`. At (or below) the critical entropy the dew
/// line meets the critical point, so `p_crit` is returned directly. Above the
/// triple-point saturated-vapour entropy it clamps to `p_min`.
#[inline]
pub fn dew_point_pressure_from_entropy(s0: SpecificHeatCapacity) -> Pressure {
    let p_min = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let p_crit = p_crit_water();

    // at / below the critical entropy the dew line terminates at the
    // critical point — return the critical pressure directly
    if s0 <= s_crit_water() { return p_crit; }

    // saturated-vapour entropy at pressure p
    let s_g = |p: Pressure| -> SpecificHeatCapacity {
        s_tp_eqm_two_phase(sat_temp_4(p), p, 1.0)
    };

    // above the triple-point saturated-vapour entropy, clamp to p_min
    if s0 >= s_g(p_min) { return p_min; }

    // seed a tight bracket from the saturation lookup table, then refine by
    // bisection within it (s_g is monotonically decreasing in p)
    let (mut p_lo, mut p_hi) = dew_point_bracket(s0);
    for _ in 0..40 {
        let p_mid = 0.5 * (p_lo + p_hi);
        if ((p_hi - p_lo) / p_mid).get::<ratio>() < 1e-9 { break; }
        if s_g(p_mid) > s0 {
            // s_g still above s0 -> need higher pressure to bring it down
            p_lo = p_mid;
        } else {
            p_hi = p_mid;
        }
    }
    0.5 * (p_lo + p_hi)
}

/// Same as [`dew_point_pressure_from_entropy`] but takes a `(p, h)` stagnation
/// state and uses its entropy. Convenient for reading superheated-vapour /
/// supercritical points straight off a p-h diagram.
#[inline]
pub fn dew_point_pressure_ph(p0: Pressure, h0: AvailableEnergy) -> Pressure {
    let s0 = s_ph_eqm(p0, h0);
    dew_point_pressure_from_entropy(s0)
}
