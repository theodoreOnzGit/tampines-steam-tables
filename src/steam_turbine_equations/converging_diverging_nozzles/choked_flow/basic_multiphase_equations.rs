use uom::si::f64::*;
use crate::interfaces::functional_programming::hs_flash_eqm::p_hs_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::mass_flux_ps_eqm_throat;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
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

