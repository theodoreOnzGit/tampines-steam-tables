use uom::ConstZero;
// note: vibe coded choked flow solution
use uom::si::f64::*;
use uom::si::pressure::bar;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::TampinesSteamTableCV;
use crate::prelude::functional_programming::ph_flash_eqm::w_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::{h_ps_eqm, v_ps_eqm};

/// Result type for choked nozzle solver
#[derive(Debug, Clone)]
pub struct ChokedNozzleSolution {
    /// Whether the flow is choked (sonic at throat)
    pub is_choked: bool,
    
    /// Actual mass flow rate (may be clamped to mdot_max if choked)
    pub mass_flowrate: MassRate,
    
    /// Maximum possible mass flow (choke limit)
    pub mdot_max: MassRate,
    
    /// Throat state (where M=1 if choked)
    pub throat_pressure: Pressure,
    pub throat_enthalpy: AvailableEnergy,
    pub throat_density: MassDensity,
    pub throat_velocity: Velocity,
    pub throat_mach: Ratio,
    
    /// Exit state
    pub exit_pressure: Pressure,
    pub exit_enthalpy: AvailableEnergy,
    pub exit_density: MassDensity,
    pub exit_velocity: Velocity,
    pub exit_mach: Ratio,
}

/// Choked nozzle solver
/// 
/// Given:
/// - Inlet stagnation-like state (p1, h1) with small v1
/// - Throat area A_star
/// - Exit area A_exit
/// - Requested mass flow OR back pressure
/// 
/// Returns: complete solution with choking status
#[inline]
pub fn solve_choked_nozzle_isentropic(
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity,
    a_star: Area,      // throat area
    a_exit: Area,      // exit area
    requested_mdot: Option<MassRate>,  // if Some: solve for exit p; if None: need p_back
    p_back: Option<Pressure>,          // back pressure (used if requested_mdot is None)
) -> Result<ChokedNozzleSolution, String> {

    // -------------------------------
    // 1. Inlet state
    // -------------------------------
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    
    let rho1 = state_1.get_specific_volume().recip();
    let s1 = state_1.get_specific_entropy();
    let s_isentrope = s1; // constant along isentrope
    
    // Stagnation enthalpy (accounting for inlet KE)
    let h0 = h1 + 0.5 * v1 * v1;
    
    // -------------------------------
    // 2. Find choke limit: max(rho*a) along isentrope at throat
    // -------------------------------
    let p_min = Pressure::new::<bar>(0.001);
    let n_scan = 500;
    
    let mut mdot_max = MassRate::ZERO;
    let mut p_star_choke = p1;
    let mut h_star_choke = h1;
    let mut rho_star_choke = rho1;
    let mut speed_of_sound_choke = Velocity::ZERO;
    
    for i in 0..=n_scan {
        // this is scanning for pressures
        // using the p,s algorithm
        let frac = (i as f64) / (n_scan as f64);
        let p = p_min + frac * (p1 - p_min);
        
        let h = h_ps_eqm(p, s_isentrope);
        let rho = v_ps_eqm(p, s_isentrope).recip();
        let speed_of_sound = w_ph_eqm(p, h); // speed of sound
        
        let mdot_candidate = rho * a_star * speed_of_sound;
        
        if mdot_candidate > mdot_max {
            mdot_max = mdot_candidate;
            p_star_choke = p;
            h_star_choke = h;
            rho_star_choke = rho;
            speed_of_sound_choke = speed_of_sound;
        }
    }
    
    // Throat velocity at choke = speed of sound
    let v_star_choke = speed_of_sound_choke;
    
    // -------------------------------
    // 3. Determine operating mode
    // -------------------------------
    let (mdot_actual, is_choked, solve_mode): (MassRate, bool, &str) = match (requested_mdot, p_back) {
        (Some(mdot_req), _) => {
            // Mode A: user specifies mdot, we find exit pressure
            if mdot_req > mdot_max {
                // Cannot achieve requested mdot isentropically
                // Clamp to choke limit
                (mdot_max, true, "mdot-specified-choked")
            } else {
                (mdot_req, false, "mdot-specified-unchoked")
            }
        }
        (None, Some(pb)) => {
            // Mode B: user specifies back pressure, we find mdot
            // Try unchoked first by computing what mdot would be needed
            // This is complex; for now we compute mdot at exit assuming isentropic to pb
            
            let h_exit_pb = h_ps_eqm(pb, s_isentrope);
            let rho_exit_pb = v_ps_eqm(pb, s_isentrope).recip();
            let v_exit_pb = (2.0 * (h0 - h_exit_pb)).sqrt();
            let mdot_pb = rho_exit_pb * a_exit * v_exit_pb;
            
            if mdot_pb > mdot_max {
                // Choked: cannot pull more than mdot_max
                (mdot_max, true, "pback-specified-choked")
            } else {
                (mdot_pb, false, "pback-specified-unchoked")
            }
        }
        (None, None) => {
            return Err("Must specify either requested_mdot or p_back".to_string());
        }
    };
    
    // -------------------------------
    // 4. Solve for exit state
    // -------------------------------
    let (p_exit, h_exit, rho_exit, v_exit): (Pressure, AvailableEnergy, MassDensity, Velocity);
    
    if is_choked {
        // Throat is sonic; downstream depends on area ratio and back pressure
        // For simplicity: if exit area > throat area (diverging section),
        // we'd need to handle supersonic expansion or shocks.
        // Here we'll do a simplified approach:
        
        if a_exit <= a_star {
            // Converging nozzle: exit is at throat (sonic)
            p_exit = p_star_choke;
            h_exit = h_star_choke;
            rho_exit = rho_star_choke;
            v_exit = v_star_choke;
        } else {
            // Converging-diverging: need to solve supersonic branch or check for shock
            // Simplified: assume isentropic expansion to exit
            // We solve: mdot = rho(p,s)*A_exit*v(p) with v from energy
            
            // Use bisection on exit pressure
            let p_exit_result = solve_exit_pressure_choked(
                h0, s_isentrope, mdot_actual, a_exit, p_star_choke, p_min
            )?;
            
            p_exit = p_exit_result;
            h_exit = h_ps_eqm(p_exit, s_isentrope);
            rho_exit = v_ps_eqm(p_exit, s_isentrope).recip();
            v_exit = (2.0 * (h0 - h_exit)).sqrt();
        }
    } else {
        // Unchoked: solve for exit state matching mdot and continuity
        
        match solve_mode {
            "mdot-specified-unchoked" => {
                // Need to find exit pressure such that continuity holds
                p_exit = solve_exit_pressure_unchoked(
                    h0, s_isentrope, mdot_actual, a_exit, p1, p_min
                )?;
                
                h_exit = h_ps_eqm(p_exit, s_isentrope);
                rho_exit = v_ps_eqm(p_exit, s_isentrope).recip();
                v_exit = (2.0 * (h0 - h_exit)).sqrt();
            }
            "pback-specified-unchoked" => {
                // Exit pressure is the back pressure
                p_exit = p_back.unwrap();
                h_exit = h_ps_eqm(p_exit, s_isentrope);
                rho_exit = v_ps_eqm(p_exit, s_isentrope).recip();
                v_exit = (2.0 * (h0 - h_exit)).sqrt();
            }
            _ => unreachable!()
        }
    }
    
    // -------------------------------
    // 5. Compute Mach numbers
    // -------------------------------
    let a_exit_sound = w_ph_eqm(p_exit, h_exit);
    let mach_exit = v_exit / a_exit_sound;
    
    let mach_throat = if is_choked {
        Ratio::new::<ratio>(1.0)
    } else {
        // Throat velocity from continuity
        let v_throat_unchoked = mdot_actual / (rho_star_choke * a_star);
        v_throat_unchoked / speed_of_sound_choke
    };
    
    // -------------------------------
    // 6. Package result
    // -------------------------------
    Ok(ChokedNozzleSolution {
        is_choked,
        mass_flowrate: mdot_actual,
        mdot_max,
        throat_pressure: p_star_choke,
        throat_enthalpy: h_star_choke,
        throat_density: rho_star_choke,
        throat_velocity: if is_choked { v_star_choke } else { mdot_actual / (rho_star_choke * a_star) },
        throat_mach: mach_throat,
        exit_pressure: p_exit,
        exit_enthalpy: h_exit,
        exit_density: rho_exit,
        exit_velocity: v_exit,
        exit_mach: mach_exit,
    })
}

/// Helper: solve for exit pressure in unchoked case
/// Finds p_exit such that: mdot = rho(p,s)*A_exit*sqrt(2*(h0-h(p,s)))
fn solve_exit_pressure_unchoked(
    h0: AvailableEnergy,
    s: SpecificHeatCapacity,
    mdot: MassRate,
    a_exit: Area,
    p_high: Pressure,
    p_low: Pressure,
) -> Result<Pressure, String> {
    let tol = 1e-6;
    let max_iter = 100;
    
    let f = |p: Pressure| -> MassRate {
        let h = h_ps_eqm(p, s);
        let rho = v_ps_eqm(p, s).recip();
        let v = (2.0 * (h0 - h)).sqrt();
        rho * a_exit * v - mdot
    };
    
    let mut p_lo = p_low;
    let mut p_hi = p_high;
    let mut f_lo = f(p_lo);
    let mut f_hi = f(p_hi);
    
    // Check bracket
    if (f_lo > MassRate::ZERO) == (f_hi > MassRate::ZERO) {
        return Err("Cannot bracket exit pressure in unchoked solver".to_string());
    }
    
    for _ in 0..max_iter {
        let p_mid = 0.5 * (p_lo + p_hi);
        let f_mid = f(p_mid);
        
        if (f_mid / mdot).get::<ratio>().abs() < tol {
            return Ok(p_mid);
        }
        
        if (f_lo > MassRate::ZERO) == (f_mid > MassRate::ZERO) {
            p_lo = p_mid;
            f_lo = f_mid;
        } else {
            p_hi = p_mid;
            f_hi = f_mid;
        }
    }
    
    Ok(0.5 * (p_lo + p_hi))
}

/// Solve unchoked nozzle: finds p_exit such that continuity, energy, AND momentum all close.
/// 
/// The three equations are:
/// 1. Continuity: mdot = rho(p2,s)*A2*v2
/// 2. Energy:     v2 = sqrt(2*(h0 - h(p2,s)))
/// 3. Momentum:   p1*A1 + mdot*v1 = p2*A2 + mdot*v2
///
/// We solve by finding p2 such that the momentum residual is zero,
/// with v2 coming from energy and rho2 from (p2,s).
fn solve_exit_pressure_unchoked_full(
    p1: Pressure,
    h0: AvailableEnergy,
    v1: Velocity,
    s: SpecificHeatCapacity,
    mdot: MassRate,
    a1: Area,
    a_exit: Area,
    p_high: Pressure,
    p_low: Pressure,
) -> Result<Pressure, String> {
    let tol = 1e-6;
    let max_iter = 100;
    
    // Momentum residual as function of p2:
    // f(p2) = (p1*A1 + mdot*v1) - (p2*A2 + mdot*v2(p2))
    // where v2(p2) comes from energy balance
    let f = |p2: Pressure| -> Force {
        let h2 = h_ps_eqm(p2, s);
        let v2_energy = (2.0 * (h0 - h2)).sqrt();
        
        let lhs = p1 * a1 + mdot * v1;
        let rhs = p2 * a_exit + mdot * v2_energy;
        
        lhs - rhs
    };
    
    let mut p_lo = p_low;
    let mut p_hi = p_high;
    let mut f_lo = f(p_lo);
    let mut f_hi = f(p_hi);
    
    // Check bracket
    if (f_lo > Force::ZERO) == (f_hi > Force::ZERO) {
        return Err("Cannot bracket exit pressure (momentum balance has no root)".to_string());
    }
    
    // Bisection
    for _ in 0..max_iter {
        let p_mid = 0.5 * (p_lo + p_hi);
        let f_mid = f(p_mid);
        
        if (f_mid / (p1 * a1)).get::<ratio>().abs() < tol {
            return Ok(p_mid);
        }
        
        if (f_lo > Force::ZERO) == (f_mid > Force::ZERO) {
            p_lo = p_mid;
            f_lo = f_mid;
        } else {
            p_hi = p_mid;
            f_hi = f_mid;
        }
    }
    
    Ok(0.5 * (p_lo + p_hi))
}

/// Choked nozzle solver with FULL momentum balance
#[inline]
pub fn solve_choked_nozzle_isentropic_full(
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity,
    a1: Area,         // inlet area
    a_star: Area,     // throat area
    a_exit: Area,     // exit area
    requested_mdot: Option<MassRate>,
    p_back: Option<Pressure>,
) -> Result<ChokedNozzleSolution, String> {

    // -------------------------------
    // 1. Inlet state
    // -------------------------------
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    
    let rho1 = state_1.get_specific_volume().recip();
    let s1 = state_1.get_specific_entropy();
    let s_isentrope = s1;
    
    let h0 = h1 + 0.5 * v1 * v1;
    
    // -------------------------------
    // 2. Find choke limit
    // -------------------------------
    let p_min = Pressure::new::<bar>(0.001);
    let n_scan = 500;
    
    let mut mdot_max = MassRate::ZERO;
    let mut p_star_choke = p1;
    let mut h_star_choke = h1;
    let mut rho_star_choke = rho1;
    let mut a_star_sound = Velocity::ZERO;
    
    for i in 0..=n_scan {
        let frac = (i as f64) / (n_scan as f64);
        let p = p_min + frac * (p1 - p_min);
        
        let h = h_ps_eqm(p, s_isentrope);
        let rho = v_ps_eqm(p, s_isentrope).recip();
        let a = w_ph_eqm(p, h);
        
        let mdot_candidate = rho * a_star * a;
        
        if mdot_candidate > mdot_max {
            mdot_max = mdot_candidate;
            p_star_choke = p;
            h_star_choke = h;
            rho_star_choke = rho;
            a_star_sound = a;
        }
    }
    
    let v_star_choke = a_star_sound;
    
    // -------------------------------
    // 3. Determine operating mode
    // -------------------------------
    let (mdot_actual, is_choked): (MassRate, bool) = match (requested_mdot, p_back) {
        (Some(mdot_req), _) => {
            if mdot_req > mdot_max {
                (mdot_max, true)
            } else {
                (mdot_req, false)
            }
        }
        (None, Some(pb)) => {
            // Compute what mdot would be needed for isentropic expansion to pb
            // using momentum balance
            
            let h_exit_pb = h_ps_eqm(pb, s_isentrope);
            let v_exit_pb = (2.0 * (h0 - h_exit_pb)).sqrt();
            
            // From momentum: p1*A1 + mdot*v1 = pb*A_exit + mdot*v_exit
            // Solve for mdot:
            let mdot_pb = (p1 * a1 - pb * a_exit) / (v_exit_pb - v1);
            
            if mdot_pb > mdot_max {
                (mdot_max, true)
            } else {
                (mdot_pb, false)
            }
        }
        (None, None) => {
            return Err("Must specify either requested_mdot or p_back".to_string());
        }
    };
    
    // -------------------------------
    // 4. Solve for exit state with MOMENTUM balance
    // -------------------------------
    let (p_exit, h_exit, rho_exit, v_exit): (Pressure, AvailableEnergy, MassDensity, Velocity);
    
    if is_choked {
        if a_exit <= a_star {
            // Converging only: exit is at throat
            p_exit = p_star_choke;
            h_exit = h_star_choke;
            rho_exit = rho_star_choke;
            v_exit = v_star_choke;
        } else {
            // C-D nozzle: solve momentum balance in supersonic region
            p_exit = solve_exit_pressure_unchoked_full(
                p1, h0, v1, s_isentrope, mdot_actual, a1, a_exit,
                p_star_choke, p_min
            )?;
            
            h_exit = h_ps_eqm(p_exit, s_isentrope);
            rho_exit = v_ps_eqm(p_exit, s_isentrope).recip();
            v_exit = (2.0 * (h0 - h_exit)).sqrt();
        }
    } else {
        // Unchoked: solve momentum balance
        p_exit = solve_exit_pressure_unchoked_full(
            p1, h0, v1, s_isentrope, mdot_actual, a1, a_exit,
            p1, p_min
        )?;
        
        h_exit = h_ps_eqm(p_exit, s_isentrope);
        rho_exit = v_ps_eqm(p_exit, s_isentrope).recip();
        v_exit = (2.0 * (h0 - h_exit)).sqrt();
    }
    
    // Verify continuity as a check
    let mdot_check = rho_exit * a_exit * v_exit;
    let continuity_error = ((mdot_check - mdot_actual) / mdot_actual).get::<ratio>().abs();
    
    if continuity_error > 0.01 {
        return Err(format!(
            "Continuity error too large: {:.2}%. Likely inconsistent inputs or need non-isentropic model.",
            continuity_error * 100.0
        ));
    }
    
    // -------------------------------
    // 5. Compute Mach numbers
    // -------------------------------
    let a_exit_sound = w_ph_eqm(p_exit, h_exit);
    let mach_exit = v_exit / a_exit_sound;
    
    let mach_throat = if is_choked {
        Ratio::new::<ratio>(1.0)
    } else {
        let v_throat = mdot_actual / (rho_star_choke * a_star);
        v_throat / a_star_sound
    };
    
    Ok(ChokedNozzleSolution {
        is_choked,
        mass_flowrate: mdot_actual,
        mdot_max,
        throat_pressure: p_star_choke,
        throat_enthalpy: h_star_choke,
        throat_density: rho_star_choke,
        throat_velocity: if is_choked { v_star_choke } else { mdot_actual / (rho_star_choke * a_star) },
        throat_mach: mach_throat,
        exit_pressure: p_exit,
        exit_enthalpy: h_exit,
        exit_density: rho_exit,
        exit_velocity: v_exit,
        exit_mach: mach_exit,
    })
}

/// Helper: solve for exit pressure in choked C-D nozzle
/// Similar to unchoked but searches in supersonic region (p < p_throat)
fn solve_exit_pressure_choked(
    h0: AvailableEnergy,
    s: SpecificHeatCapacity,
    mdot: MassRate,
    a_exit: Area,
    p_throat: Pressure,
    p_min: Pressure,
) -> Result<Pressure, String> {
    // In supersonic region: p_exit < p_throat
    solve_exit_pressure_unchoked(h0, s, mdot, a_exit, p_throat, p_min)
}



