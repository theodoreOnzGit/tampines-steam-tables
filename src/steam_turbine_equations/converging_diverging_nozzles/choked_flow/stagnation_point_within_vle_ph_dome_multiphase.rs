// Critical-flow solvers for the case where the STAGNATION state lies
// inside the p-h VLE dome (two-phase, at or below the critical point).
//
// These exploit the simplification that an isentrope starting inside the
// dome stays inside it on depressurisation (the dome only widens as p
// falls), so no region switching or flashing event needs handling here.

use uom::ConstZero;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::pressure::megapascal;
use uom::si::pressure::pascal;

use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;

/// Critical pressure & mass flux for a stagnation state that sits
/// INSIDE the p-h VLE dome (two-phase, at or below the critical point).
///
/// Precondition: (p0, h0) is two-phase — i.e. ph_flash_region(p0,h0) == Region4.
/// Once inside the dome, isentropic depressurisation stays inside it
/// (the dome only widens as p falls), so there is no flashing event and
/// no region switching to handle here.
///
/// Method (Moody / max-flux form of the HEM choking criterion):
///   along the isentrope s = s0,
///     G(p) = rho(p,s0) * sqrt( 2 * (h0 - h(p,s0)) )
///   G(p0) = 0, rises to a single interior maximum at the choke point,
///   then falls as rho -> 0. The choke is argmax_p G(p).
///
/// This avoids mass_flux_ps_eqm_throat (finite-difference sound speed +
/// bubble-point clamp) entirely; it only needs smooth h(p,s0), v(p,s0).
///
/// Consistent with the validated inverse map: max-G <=> Mach 1 <=>
/// h0 = h_t + 0.5 * u_t^2  (get_stagnation_conditions_from_throat_ps).
#[inline]
pub fn get_critical_pressure_and_mass_flux_ph_vle_dome(
    p0: Pressure,
    h0: AvailableEnergy,
) -> (Pressure, MassFlux) {

    // isentrope to march down
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

    // Maximise G over [p_min, p0] by golden-section search. G is unimodal
    // (zero at p0, single interior peak at the choke, falling toward p_min),
    // so this is robust and needs no derivative of the noisy sound speed.
    //
    // gr = (sqrt(5) - 1)/2 ~= 0.618 is the golden ratio: it places the two
    // interior probes so that one probe is reused after each bracket
    // reduction, costing one G-evaluation per iteration.
    //
    // Reference:
    //   Price, C. J., & Robertson, B. L. (2012). Golden Section Search.
    //   In Encyclopedia of Engineering Optimization and Heuristics
    //   (pp. 1-4). Singapore: Springer Nature Singapore.
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;     // 0.618...
    let mut a = p_min.get::<pascal>();
    let mut b = p0.get::<pascal>();
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);

    let max_iter = 100;
    let tol_pa = 1.0;                          // 1 Pa bracket width is plenty
    for _ in 0..max_iter {
        if (b - a).abs() < tol_pa { break; }
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

    let p_crit = Pressure::new::<pascal>(0.5 * (a + b));
    let g_crit = g_of_p(p_crit.get::<pascal>());

    (p_crit, g_crit)
}
