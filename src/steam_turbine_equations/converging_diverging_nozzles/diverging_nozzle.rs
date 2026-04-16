use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::pressure::pascal;
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
    s0: SpecificHeatCapacity,
    p2: Pressure,
    a_throat: Area,
    a_exit: Area,
    mass_rate_throat: MassRate,
    state_throat: TampinesSteamTableCV,
) -> (Velocity, TampinesSteamTableCV) {
    
    // Calculate reference mass flux (must be conserved through nozzle)
    let mass_flux_ref: MassFlux = mass_rate_throat / a_throat;
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let inlet_stagnation_state = 
        TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);

    // Calculate perfectly expanded solution 
    let (p_ideal_expansion, v_ideal_expansion, state_ideal_expansion) 
        = calculate_isentropic_exit_pressure_velocity_and_state(
            inlet_stagnation_state, 
            a_exit, 
            mass_rate_throat,
        );

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

    // firstly, we want to check the nozzle boundary pressure,
    // it should be at least, the ideal expansion pressure
    // if p2 is lower than this ideal expansion pressure, oblique shocks 
    // will form outside
    let mut p2_nozzle_boundary = p_ideal_expansion;
    
    // For isentropic flow: s₂ = s_throat
    let s2_isentropic = state_throat.get_specific_entropy();
    
    // (p,s) flash to get isentropic outlet state
    let state_2_isentropic = 
        TampinesSteamTableCV::new_from_ps(
            p2_nozzle_boundary, s2_isentropic, ref_vol
        );
    let h2_isentropic = state_2_isentropic.get_specific_enthalpy();
    
    // Check if isentropic solution satisfies mass balance
    let mass_flux_isentropic = 
        calculate_mass_flux_at_outlet(h0, p2_nozzle_boundary, h2_isentropic);
    
    let mass_flux_error: f64 = 
        ((mass_flux_isentropic - mass_flux_ref) / mass_flux_ref).get::<ratio>();
    
    const TOLERANCE: f64 = 0.0001;  // 0.01% tolerance
    
    let pressure_tolerance = Pressure::new::<pascal>(100.0); // or appropriate tolerance
    let pressure_diff = (p2 - p_ideal_expansion).abs();

    if mass_flux_error.abs() < TOLERANCE && pressure_diff < pressure_tolerance {
        // Isentropic solution is valid!
        // That means either we have perfect expansions
        let v_outlet: Velocity = v_ideal_expansion;
        let state_outlet = state_ideal_expansion;
        
        return (v_outlet, state_outlet);
    }
    if mass_flux_error.abs() < TOLERANCE && p2 < p_ideal_expansion {
        // if outlet pressure is more than ideal expansion pressure, we 
        // will have the correct mass flux in the outlet
        // in this case, we will have oblique shocks outside the nozzle
        let h_nozzle_outlet = h2_isentropic;
        let v_nozzle_outlet: Velocity = v_ideal_expansion;
        //let state_nozzle_outlet = state_ideal_expansion;

        // now after this ideal expansion, 
        // we should have a certain enthalpy and entropy
        //
        // (p_nozzle_outlet, s_ideal) -> (p2, unknown state)
        //
        // Note: (p_nozzle_outlet > p2)
        //
        // I'm not quite sure as to how expansion is going to occur
        // But there is going to be further pressure decrease, and then 
        // some mixing
        //
        // What is going to be our state after these shocks?
        // Indeed, in the most ideal case, it is further isentropic 
        // expansion to achieve higher velocities 
        //
        // In the non ideal case, we assume there is not velocity increase, 
        // but a pressure decrease, ie joule thompson effect. 
        // that is after mixing and such
        //
        let state_outlet = TampinesSteamTableCV::new_from_ph(
            p2, h_nozzle_outlet, ref_vol
        );
        // this will give some estimate as to what the outlet state should 
        // be. A conservative estimate
        //
        // We won't be doing a mass conservation equation so to speak.


        return (v_nozzle_outlet, state_outlet);
    }

    // ========================================================================
    // Step 2: Non-isentropic solution (shocks present) - Use bisection
    // ========================================================================
    
    // Physical bounds on outlet enthalpy:
    // - Lower bound: h2_isentropic (minimum possible, maximum expansion)
    // - Upper bound: h0 (maximum possible, zero velocity)
    let mut h_lower = h2_isentropic;
    let mut h_upper = h0;
    p2_nozzle_boundary = p2;
    
    let max_iterations = 50;
    let enthalpy_tolerance = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0);
    
    // Bisection loop to find h₂ that satisfies mass balance
    for _iteration in 0..max_iterations {
        // Midpoint guess
        let h_mid = 0.5 * (h_lower + h_upper);
        
        // Calculate mass flux at this enthalpy
        let mass_flux_guess = calculate_mass_flux_at_outlet(h0, p2, h_mid);
        
        // Check error
        let error: f64 = 
            ((mass_flux_guess - mass_flux_ref) / mass_flux_ref).get::<ratio>();
        dbg!(&mass_flux_guess);
        dbg!(&mass_flux_ref);
        dbg!(&error);
        
        // Check if converged
        if error.abs() < TOLERANCE {
            let h_outlet = h_mid;
            let v_outlet: Velocity = (2.0 * (h0 - h_outlet)).sqrt();
            let state_outlet = TampinesSteamTableCV::new_from_ph(
                p2_nozzle_boundary, h_outlet, ref_vol
            );
            
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
            let state_outlet = TampinesSteamTableCV::new_from_ph(
                p2_nozzle_boundary, h_outlet, ref_vol
            );
            dbg!(&(mass_flux_ref,mass_flux_guess));
            dbg!(&(h_upper,h_lower));
            
            return (v_outlet, state_outlet);
        }
    }
    
    // ========================================================================
    // Step 3: Max iterations reached - return best guess
    // ========================================================================
    
    let h_outlet = 0.5 * (h_lower + h_upper);
    let v_outlet: Velocity = (2.0 * (h0 - h_outlet)).sqrt();
    let state_outlet = TampinesSteamTableCV::new_from_ph(
        p2_nozzle_boundary, h_outlet, ref_vol
    );
    
    return (v_outlet, state_outlet);
}

/// Calculate exit pressure for isentropic expansion through CD nozzle
/// assuming choked flow
///
/// this is for perfectly expanded flow
#[inline]
pub fn calculate_isentropic_exit_pressure_velocity_and_state(
    inlet_stagnation_state: TampinesSteamTableCV,
    a_exit: Area,
    mass_flowrate_choked: MassRate,
) -> (Pressure, Velocity, TampinesSteamTableCV) {
    
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    
    // For isentropic flow: s_exit = s0
    // Mass continuity: ṁ = ρ_exit * v_exit * A_exit
    // Energy: v_exit = sqrt(2*(h0 - h_exit))
    //
    // Need to find p_exit such that these are satisfied
    let h0: AvailableEnergy = inlet_stagnation_state.get_specific_enthalpy();
    let s0: SpecificHeatCapacity = inlet_stagnation_state.get_specific_entropy();
    let p0: Pressure = inlet_stagnation_state.get_pressure();
    
    // Use bisection to find p_exit
    let mut p_lower = Pressure::new::<pascal>(1000.0);  // Very low pressure
    let mut p_upper = p0;      
    
    let max_iterations = 50;
    let tolerance = Pressure::new::<pascal>(100.0);
    let mut state_exit: TampinesSteamTableCV;
    let mut v_exit: Velocity;
    
    for _ in 0..max_iterations {
        let p_mid = 0.5 * (p_lower + p_upper);
        
        // Calculate state at this pressure (isentropic)
        state_exit = TampinesSteamTableCV::new_from_ps(p_mid, s0, ref_vol);
        let h_exit = state_exit.get_specific_enthalpy();
        let rho_exit = state_exit.get_rho();
        
        // Calculate velocity from energy equation
        v_exit = (2.0 * (h0 - h_exit)).sqrt();
        
        // Calculate mass flowrate
        let mass_flowrate_calc = rho_exit * v_exit * a_exit;
        
        // Check error
        let error = (mass_flowrate_calc - mass_flowrate_choked) / mass_flowrate_choked;
        
        if error.get::<ratio>().abs() < 0.0001 {
            return (p_mid, v_exit, state_exit);
        }
        
        // Adjust bounds
        // Lower pressure → higher velocity → higher mass flow (for supersonic)
        if error.get::<ratio>() > 0.0 {
            // Mass flow too high, increase pressure
            p_lower = p_mid;
        } else {
            // Mass flow too low, decrease pressure
            p_upper = p_mid;
        }
        
        if (p_upper - p_lower) < tolerance {
            return (p_mid, v_exit, state_exit);
        }
    }

    let p_mid = 0.5 * (p_lower + p_upper);
    state_exit = TampinesSteamTableCV::new_from_ps(p_mid, s0, ref_vol);
    let h_exit = state_exit.get_specific_enthalpy();
    v_exit = (2.0 * (h0 - h_exit)).sqrt();
    
    return (p_mid, v_exit, state_exit);
}
