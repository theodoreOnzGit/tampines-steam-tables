use uom::ConstZero;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::pressure::megapascal;
use uom::si::pressure::pascal;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;

use crate::constants::p_crit_water;
use crate::constants::t_crit_water;
use crate::interfaces::functional_programming::hs_flash_eqm::p_hs_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::t_ph_eqm;
use crate::interfaces::functional_programming::ps_flash_eqm::ps_flash_region;
use crate::interfaces::functional_programming::ps_flash_eqm::w_ps_wood_wallis;
use crate::interfaces::functional_programming::ps_flash_eqm::x_ps_flash;
use crate::prelude::functional_programming::ps_flash_eqm::mass_flux_ps_eqm_throat;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::prelude::functional_programming::pt_flash_eqm::FwdEqnRegion;

use crate::region_1_subcooled_liquid::s_tp_1;
use crate::region_1_subcooled_liquid::v_tp_1;
use crate::region_2_vapour::s_tp_2;
use crate::region_2_vapour::v_tp_2;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
use crate::steam_turbine_equations::choked_flow::single_phase_basic_choked_flow::get_critical_pressure_pure_vapour_ph_stagnation_properties;

/// these contain choked flow algorithms for single phase choked flow,
///
/// whether be it finding critical pressure for ideal gas, or for those 
/// where the choked flow is in the pure vapour phase
pub mod single_phase_basic_choked_flow;



/// gets critical pressure and mass flux for water and steam 
/// given stagnation properties,
/// should work for all given regions of steam table
///
/// Note, this relies on the ph, ps algorithm, for region 5, it is 
/// not thoroughly implemented yet
/// gets critical pressure and mass flux for water and steam 
/// given stagnation properties,
/// should work for all given regions of steam table
///
/// Note, this relies on the ph, ps algorithm, for region 5, it is 
/// not thoroughly implemented yet
/// gets critical pressure and mass flux for water and steam 
/// given stagnation properties,
/// should work for all given regions of steam table
///
/// Note, this relies on the ph, ps algorithm, for region 5, it is 
/// not thoroughly implemented yet
#[inline]
pub fn get_critical_pressure_and_mass_flux_with_stagnation_props(
    s0: SpecificHeatCapacity,
    h0: AvailableEnergy,
    p0: Pressure) -> (Pressure, MassFlux) {

    let debug = true;

    // now before anything, we want to get the region of the scans 
    let region_stagnation_props = ph_flash_region(p0, h0);

    // for high entropy states (s0 >= 9.2 kJ/kg/K), 
    // isentropic depressurisation stays in single phase vapour
    // so we can use the pure vapour algorithm directly
    if s0 >= SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.2) {
        if debug { println!("using vapour only algorithm, s0 >= 9.2 kJ/kg/K"); }

        let s0_opt = Some(s0);
        let critical_pressure_choked_flow 
            = get_critical_pressure_pure_vapour_ph_stagnation_properties(
                p0, h0, s0_opt
            );

        let c = w_ps_wood_wallis(critical_pressure_choked_flow, s0);
        let rho_throat = v_ps_eqm(critical_pressure_choked_flow, s0).recip();
        let critical_mass_flux = c * rho_throat;

        if debug {
            dbg!(&(s0, critical_pressure_choked_flow, critical_mass_flux));
        }

        return (critical_pressure_choked_flow, critical_mass_flux);
    }

    // for region 2 superheated vapour with s0 < 9.2 kJ/kg/K,
    // check if isentrope crosses into Region 4 during depressurisation
    // by comparing s0 with s_g at stagnation pressure
    //
    // if s0 >= s_g(p0), isentrope stays in vapour — use pure vapour algorithm
    // if s0 < s_g(p0), isentrope crosses into Region 4 — use generalised algorithm
    if region_stagnation_props == FwdEqnRegion::Region2 {
        let p_crit_water_val = p_crit_water();
        if p0 < p_crit_water_val {
            let t_sat = sat_temp_4(p0);
            let s_g = s_tp_2(t_sat, p0);
            if s0 >= s_g {
                if debug { println!("Region 2, s0 >= s_g, using pure vapour algorithm"); }
                let s0_opt = Some(s0);
                let critical_pressure_choked_flow 
                    = get_critical_pressure_pure_vapour_ph_stagnation_properties(
                        p0, h0, s0_opt
                    );
                let c = w_ps_wood_wallis(critical_pressure_choked_flow, s0);
                let rho_throat = v_ps_eqm(critical_pressure_choked_flow, s0).recip();
                let critical_mass_flux = c * rho_throat;
                if debug {
                    dbg!(&(s0, s_g, critical_pressure_choked_flow, critical_mass_flux));
                }
                return (critical_pressure_choked_flow, critical_mass_flux);
            }
            if debug { println!("Region 2, s0 < s_g, using generalised algorithm"); }
        }
    }

    // for region 3, use vapour algorithm only if above critical temperature
    if region_stagnation_props == FwdEqnRegion::Region3 {
        let t_crit = t_crit_water();
        let t0 = t_ph_eqm(p0, h0);
        if t0 > t_crit {
            if debug { println!("Region 3, T > T_crit, using pure vapour algorithm"); }
            let s0_opt = Some(s0);
            let critical_pressure_choked_flow 
                = get_critical_pressure_pure_vapour_ph_stagnation_properties(
                    p0, h0, s0_opt
                );
            let c = w_ps_wood_wallis(critical_pressure_choked_flow, s0);
            let rho_throat = v_ps_eqm(critical_pressure_choked_flow, s0).recip();
            let critical_mass_flux = c * rho_throat;
            if debug {
                dbg!(&(s0, critical_pressure_choked_flow, critical_mass_flux));
            }
            return (critical_pressure_choked_flow, critical_mass_flux);
        }
        if debug { println!("Region 3, T <= T_crit, using generalised algorithm"); }
    }

    // generalised algorithm for:
    // - Region 1 (subcooled liquid)
    // - Region 4 (wet steam)
    // - Region 2 with s0 < s_g (will flash during depressurisation)
    // - Region 3 liquid-like (T <= T_crit)
    //
    // root finding: G_energy = G_HEM
    // f(p) = G_energy(p) - G_HEM(p) = 0
    // positive when g_energy > g_hem (past choked point, too low pressure)
    // negative when g_hem > g_energy (not yet choked, too high pressure)

    let p_min_steam_table = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let p_decrement: Pressure = 0.05 * (p0 - p_min_steam_table);

    let root_fn = |p_test_pa: f64| -> f64 {
        let p_test = Pressure::new::<pascal>(p_test_pa);
        let region = ps_flash_region(p_test, s0);

        // no choked flow in subcooled liquid
        if region == FwdEqnRegion::Region1 {
            return -1e10_f64;
        }

        let h_test = h_ps_eqm(p_test, s0);
        let kinetic_energy = h0 - h_test;

        // pressure too low, expansion exceeded stagnation enthalpy
        if kinetic_energy < AvailableEnergy::ZERO {
            return 1e10_f64;
        }

        let rho = v_ps_eqm(p_test, s0).recip();
        let g_energy: MassFlux = rho * (2.0 * kinetic_energy).sqrt();
        let g_hem: MassFlux = mass_flux_ps_eqm_throat(p_test, s0);

        if debug {
            let quality = x_ps_flash(p_test, s0);
            dbg!(&(p_test, region, quality, g_energy, g_hem));
        }

        (g_energy - g_hem).get::<kilogram_per_square_meter_second>()
    };

    // scan downward from p0 to find sign change bracket
    let mut p_scan = p0;
    let mut f_prev = root_fn(p_scan.get::<pascal>());
    let mut p_upper_bracket = p0;
    let mut p_lower_bracket = p_min_steam_table;
    let mut found_bracket = false;

    while (p_scan - p_decrement) > p_min_steam_table {
        p_scan -= p_decrement;
        let f_curr = root_fn(p_scan.get::<pascal>());

        if debug {
            dbg!(&(p_scan, f_prev, f_curr));
        }

        // sign change detected — bracket found
        if f_prev * f_curr < 0.0 {
            p_upper_bracket = p_scan + p_decrement;
            p_lower_bracket = p_scan;
            found_bracket = true;
            if debug {
                println!("bracket found!");
                dbg!(&(p_upper_bracket, p_lower_bracket));
            }
            break;
        }

        f_prev = f_curr;
    }

    if !found_bracket {
        panic!("unable to find bracket for critical mass flux root finding");
    }

    // in-house regula falsi (false position) method
    // converges faster than bisection for smooth functions
    let max_iterations = 100;
    let tolerance = 1e-6_f64; // kg/(m²·s)

    let mut p_low = p_lower_bracket;
    let mut p_high = p_upper_bracket;
    let mut f_low = root_fn(p_low.get::<pascal>());
    let mut f_high = root_fn(p_high.get::<pascal>());
    let mut p_crit = p_low;

    for _ in 0..max_iterations {
        let dp = p_high - p_low;
        let df = f_high - f_low;

        // guard against division by zero
        if df.abs() < 1e-30 {
            break;
        }

        // regula falsi interpolation:
        // p_new = p_low - f_low * (p_high - p_low) / (f_high - f_low)
        p_crit = p_low - Pressure::new::<pascal>(
            f_low * dp.get::<pascal>() / df
        );

        // clamp to bounds
        if p_crit < p_low { p_crit = p_low; }
        if p_crit > p_high { p_crit = p_high; }

        let f_crit = root_fn(p_crit.get::<pascal>());

        if debug {
            dbg!(&(p_crit, f_crit));
        }

        // check convergence
        if f_crit.abs() < tolerance {
            break;
        }

        // update bracket
        if f_low * f_crit < 0.0 {
            p_high = p_crit;
            f_high = f_crit;
        } else {
            p_low = p_crit;
            f_low = f_crit;
        }
    }

    let g_crit = mass_flux_ps_eqm_throat(p_crit, s0);

    if debug {
        let region = ps_flash_region(p_crit, s0);
        let quality = x_ps_flash(p_crit, s0);
        dbg!(&(p_crit, region, quality, g_crit));
    }

    (p_crit, g_crit)
}
#[inline]
pub fn isentropic_pressure_scan_of_mass_flux(
    s0: SpecificHeatCapacity,
    p0: Pressure) -> () {

    let p_min_steam_table = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    
    // number of scan steps
    let n_steps = 1000;
    let dp = (p0 - p_min_steam_table) / n_steps as f64;

    let mut max_g_hem = MassFlux::ZERO;
    let mut p_crit = p0;

    let mut p_test = p0;

    for _ in 0..n_steps {
        p_test -= dp;

        if p_test < p_min_steam_table {
            break;
        }

        //// skip region 1 — no choked flow in subcooled liquid
        let region = ps_flash_region(p_test, s0);
        //if region == FwdEqnRegion::Region1 {
        //    continue;
        //}

        let g_hem = mass_flux_ps_eqm_throat(p_test, s0);

        let quality = x_ps_flash(p_test, s0);
        dbg!(&(p_test, region, quality, g_hem, max_g_hem, p_crit));

        if g_hem > max_g_hem {
            max_g_hem = g_hem;
            p_crit = p_test;
            dbg!(&("new maximum found!", p_crit, max_g_hem));
        }
    }

}





/// Analytical HEM critical mass flux from Saha (1978) NUREG/CR-0417 eq. 10
///
/// G²_max = -1 / (dv_mix/dP|_s)
///
/// where dv_mix/dP|_s is expanded as:
/// x * (dv_g/dP)_s + (v_g - v_f) * (dx/dP)_s + (1-x) * (dv_f/dP)_s
///
/// This is the analytical version of mass_flux_ps_eqm_throat
/// which computes the same quantity via finite difference
///
/// Takes throat conditions (p, s) or (p, h) — NOT stagnation conditions
///
/// This uses region 1 and 2 eqns
#[inline]
pub fn g_max_hem_analytical_ps(
    p: Pressure,
    s: SpecificHeatCapacity,
) -> MassFlux {

    let p_min = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let dp = p * 1e-5_f64;
    let p_plus  = p + dp;
    let p_minus = if p - dp > p_min { p - dp } else { p_min };
    let dp_actual = p_plus - p_minus;

    let t_sat       = sat_temp_4(p);
    let t_sat_plus  = sat_temp_4(p_plus);
    let t_sat_minus = sat_temp_4(p_minus);

    // quality at throat
    let x = x_ps_flash(p, s);

    // --- term 1: x * (dv_g/dP)_s ---
    // dv_g/dP along saturation curve
    let v_g_plus  = v_tp_2(t_sat_plus,  p_plus);
    let v_g_minus = v_tp_2(t_sat_minus, p_minus);
    let dv_g_dp = (v_g_plus - v_g_minus) / dp_actual;
    let term1 = x * dv_g_dp;

    // --- term 2: (v_g - v_f) * (dx/dP)_s ---
    // dx/dP along isentrope at constant s
    // x = (s - s_f) / (s_g - s_f)
    // so dx/dP = d/dP [(s - s_f) / (s_g - s_f)]
    //          = [-(ds_f/dP)(s_g - s_f) - (s - s_f)(ds_g/dP - ds_f/dP)]
    //            / (s_g - s_f)²
    let s_f       = s_tp_1(t_sat,       p);
    let s_f_plus  = s_tp_1(t_sat_plus,  p_plus);
    let s_f_minus = s_tp_1(t_sat_minus, p_minus);

    let s_g       = s_tp_2(t_sat,       p);
    let s_g_plus  = s_tp_2(t_sat_plus,  p_plus);
    let s_g_minus = s_tp_2(t_sat_minus, p_minus);

    let ds_f_dp = (s_f_plus - s_f_minus) / dp_actual;
    let ds_g_dp = (s_g_plus - s_g_minus) / dp_actual;

    let s_fg = s_g - s_f;
    let dx_dp = (-ds_f_dp * s_fg - (s - s_f) * (ds_g_dp - ds_f_dp))
               / (s_fg * s_fg);

    let v_g = v_tp_2(t_sat, p);
    let v_f = v_tp_1(t_sat, p);
    let term2 = (v_g - v_f) * dx_dp;

    // --- term 3: (1-x) * (dv_f/dP)_s ---
    // dv_f/dP along saturation curve
    let v_f_plus  = v_tp_1(t_sat_plus,  p_plus);
    let v_f_minus = v_tp_1(t_sat_minus, p_minus);
    let dv_f_dp = (v_f_plus - v_f_minus) / dp_actual;
    let term3 = (1.0 - x) * dv_f_dp;

    // --- G²_max = -1 / (term1 + term2 + term3) ---
    let dv_mix_dp = term1 + term2 + term3;

    // dv_mix_dp should be negative in two-phase region
    // G_max = sqrt(-1 / dv_mix_dp)
    let g_max_squared = dv_mix_dp.recip() * -1.0;

    g_max_squared.sqrt()
}

/// same as g_max_hem_analytical_ps but takes (p, h) as input
/// converts h to s internally
#[inline]
pub fn g_max_hem_analytical_ph(
    p: Pressure,
    h: AvailableEnergy,
) -> MassFlux {
    let s = s_ph_eqm(p, h);
    g_max_hem_analytical_ps(p, s)
}

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
