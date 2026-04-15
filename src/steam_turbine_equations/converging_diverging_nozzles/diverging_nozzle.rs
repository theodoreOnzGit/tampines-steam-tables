use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::TampinesSteamTableCV;

/// given a sonic flow, 
///
/// note, shocks may occur here 
/// 
/// given a pressure at the outlet, p2,
/// and throat state, guess the state of flow going out
/// mass flowrate is based on choked flow
///
/// stagnation properties should also be supplied to facilitate calculation
///
/// note that this is no longer isentropic
#[inline]
pub fn guess_velocity_and_state_for_diverge_nozzle_from_choked_throat(
    h0: AvailableEnergy,
    p2: Pressure,
    a_throat: Area,
    mass_rate_throat: MassRate,
    state_throat: TampinesSteamTableCV,
) -> (Velocity, TampinesSteamTableCV) {
    
    // Calculate reference mass flux (must be conserved through nozzle)
    let mass_flux_ref: MassFlux = mass_rate_throat / a_throat;
    let ref_vol = Volume::new::<cubic_meter>(1.0);

    // Helper: Calculate mass flux given outlet enthalpy (p,h) flash
    fn calculate_mass_flux_at_outlet(
        h0: AvailableEnergy,
        p2: Pressure,
        h2: AvailableEnergy,
    ) -> MassFlux {
        // Energy equation: v₂ = √(2(h₀ - h₂))
        let v2: Velocity = (2.0 * (h0 - h2)).sqrt();
        
        // Get density from (p,h) flash
        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_2 = TampinesSteamTableCV::new_from_ph(p2, h2, ref_vol);
        let rho2 = state_2.get_rho();
        
        // Mass flux: G = ρv
        let mass_flux: MassFlux = rho2 * v2;
        
        mass_flux
    }

    // ========================================================================
    // Step 1: Try isentropic solution (no shocks)
    // ========================================================================
    
    // For isentropic flow: s₂ = s_throat
    let s2_isentropic = state_throat.get_specific_entropy();
    
    // (p,s) flash to get isentropic outlet state
    let state_2_isentropic = 
        TampinesSteamTableCV::new_from_ps(p2, s2_isentropic, ref_vol);
    let h2_isentropic = state_2_isentropic.get_specific_enthalpy();
    
    // Check if isentropic solution satisfies mass balance
    let mass_flux_isentropic = 
        calculate_mass_flux_at_outlet(h0, p2, h2_isentropic);
    
    let mass_flux_error: f64 = 
        ((mass_flux_isentropic - mass_flux_ref) / mass_flux_ref).get::<ratio>();
    
    const TOLERANCE: f64 = 0.0001;  // 0.01% tolerance
    
    if mass_flux_error.abs() < TOLERANCE {
        // Isentropic solution is valid!
        let h_outlet = h2_isentropic;
        let v_outlet: Velocity = (2.0 * (h0 - h_outlet)).sqrt();
        let state_outlet = TampinesSteamTableCV::new_from_ph(p2, h_outlet, ref_vol);
        
        return (v_outlet, state_outlet);
    }

    // ========================================================================
    // Step 2: Non-isentropic solution (shocks present) - Use bisection
    // ========================================================================
    
    // Physical bounds on outlet enthalpy:
    // - Lower bound: h2_isentropic (minimum possible, maximum expansion)
    // - Upper bound: h0 (maximum possible, zero velocity)
    let mut h_lower = h2_isentropic;
    let mut h_upper = h0;
    
    let max_iterations = 50;
    let enthalpy_tolerance = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0);
    
    // Bisection loop to find h₂ that satisfies mass balance
    for iteration in 0..max_iterations {
        // Midpoint guess
        let h_mid = 0.5 * (h_lower + h_upper);
        
        // Calculate mass flux at this enthalpy
        let mass_flux_guess = calculate_mass_flux_at_outlet(h0, p2, h_mid);
        
        // Check error
        let error: f64 = 
            ((mass_flux_guess - mass_flux_ref) / mass_flux_ref).get::<ratio>();
        
        // Check if converged
        if error.abs() < TOLERANCE {
            let h_outlet = h_mid;
            let v_outlet: Velocity = (2.0 * (h0 - h_outlet)).sqrt();
            let state_outlet = TampinesSteamTableCV::new_from_ph(p2, h_outlet, ref_vol);
            
            return (v_outlet, state_outlet);
        }
        
        // Adjust bounds based on error
        // Physical reasoning: higher h₂ → lower v₂ → lower mass flux
        if error > 0.0 {
            // Mass flux too high, need to increase h₂
            h_lower = h_mid;
        } else {
            // Mass flux too low, need to decrease h₂
            h_upper = h_mid;
        }
        
        // Check if bounds have converged
        if (h_upper - h_lower) < enthalpy_tolerance {
            let h_outlet = 0.5 * (h_lower + h_upper);
            let v_outlet: Velocity = (2.0 * (h0 - h_outlet)).sqrt();
            let state_outlet = TampinesSteamTableCV::new_from_ph(p2, h_outlet, ref_vol);
            
            return (v_outlet, state_outlet);
        }
    }
    
    // ========================================================================
    // Step 3: Max iterations reached - return best guess
    // ========================================================================
    
    let h_outlet = 0.5 * (h_lower + h_upper);
    let v_outlet: Velocity = (2.0 * (h0 - h_outlet)).sqrt();
    let state_outlet = TampinesSteamTableCV::new_from_ph(p2, h_outlet, ref_vol);
    
    return (v_outlet, state_outlet);
}
