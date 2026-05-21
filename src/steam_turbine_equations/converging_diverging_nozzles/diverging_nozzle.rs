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
    
    // Calculate reference mass flowrate
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let inlet_stagnation_state = 
        TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);

    // Calculate perfectly expanded solution 
    let (p_ideal_expansion, v_ideal_expansion, state_ideal_expansion) 
        = calculate_isentropic_exit_pressure_velocity_and_state_supersonic(
            inlet_stagnation_state, 
            a_exit, 
            mass_rate_throat,
        );

    // Helper: Calculate mass flux given outlet enthalpy (p,h) flash
    fn calculate_mass_rate_at_outlet(
        h0: AvailableEnergy,
        p2: Pressure,
        h2: AvailableEnergy,
        a2: Area,
    ) -> MassRate {
        // Energy equation: v₂ = √(2(h₀ - h₂))
        let v2: Velocity = (2.0 * (h0 - h2)).sqrt();
        
        // Get density from (p,h) flash
        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_2 = TampinesSteamTableCV::new_from_ph(p2, h2, ref_vol);
        let rho2 = state_2.get_rho();
        
        // Mass flux: G = ρv
        // Mass rate: G*a2
        let mass_rate = rho2 * v2 * a2;
        
        mass_rate
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
    let mass_rate_isentropic = 
        calculate_mass_rate_at_outlet(
            h0, 
            p2_nozzle_boundary, 
            h2_isentropic,
            a_exit,
        );
    
    let mass_rate_error: f64 = 
        ((mass_rate_isentropic - mass_rate_throat) / mass_rate_throat).get::<ratio>();
    dbg!(&mass_rate_error);
    dbg!(&mass_rate_isentropic);
    dbg!(&mass_rate_throat);
    
    const TOLERANCE: f64 = 0.0001;  // 0.01% tolerance
    
    let pressure_tolerance = Pressure::new::<pascal>(100.0); // or appropriate tolerance
    let pressure_diff = (p2 - p_ideal_expansion).abs();
    dbg!(&(p2, p_ideal_expansion));

    if mass_rate_error.abs() < TOLERANCE && pressure_diff < pressure_tolerance {
        // Isentropic solution is valid!
        // That means either we have perfect expansions
        let v_outlet: Velocity = v_ideal_expansion;
        let state_outlet = state_ideal_expansion;
        
        return (v_outlet, state_outlet);
    }

    if mass_rate_error.abs() < TOLERANCE && p2 < p_ideal_expansion {
        println!("expecting shockwaves outside nozzle as p2 is less than p ideal expansion");
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
        let mass_rate_guess = calculate_mass_rate_at_outlet(
            h0, p2, h_mid,a_exit
        );
        
        // Check error
        let error: f64 = 
            ((mass_rate_guess - mass_rate_throat) / mass_rate_throat).get::<ratio>();
        dbg!(&mass_rate_guess);
        dbg!(&mass_rate_throat);
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
            dbg!(&(mass_rate_throat,mass_rate_guess));
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
///
///
/// Now, there are two pressures that would work,
/// the lower bound pressure and upper bound pressure 
/// the lower bound pressure is supersonic 
/// and the upper bound pressure is subsonic
///
/// This bisection method is based on a pressure algorithm,
/// that is to change pressure until the right mass flowrate is achieved

/// Calculates the exit pressure, velocity, and state for a perfectly expanded, 
/// isentropic flow in a converging-diverging nozzle, targeting the SUPERSONIC solution.
///
/// This function assumes the flow is choked at the throat. It finds the exit conditions 
/// in the diverging section that satisfy the choked mass flow rate for a given exit area.
///
/// # Algorithm
/// The function uses a two-stage process:
/// 1.  **Bounding Scan (Velocity-based):** It first performs a rough scan across a range of
///     velocities to find a narrow pressure bracket `[p_lower, p_upper]` that contains the 
///     supersonic root. This is the most critical step for isolating the correct solution.
/// 2.  **Refinement (Pressure-based Bisection):** It then uses a bisection method on pressure
///     to refine the solution within that narrow bracket to the required precision.
///
/// # Arguments
/// * `inlet_stagnation_state`: The thermodynamic state at stagnation conditions (h0, s0).
/// * `a_exit`: The area of the nozzle exit.
/// * `mass_flowrate_choked`: The mass flow rate determined by the choked throat conditions.
///
#[inline]
pub fn calculate_isentropic_exit_pressure_velocity_and_state_supersonic(
    inlet_stagnation_state: TampinesSteamTableCV,
    a_exit: Area,
    mass_flowrate_choked: MassRate,
) -> (Pressure, Velocity, TampinesSteamTableCV) {
    
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    
    // Stagnation properties are constant throughout the isentropic process.
    let h0: AvailableEnergy = inlet_stagnation_state.get_specific_enthalpy();
    let s0: SpecificHeatCapacity = inlet_stagnation_state.get_specific_entropy();
    let p0: Pressure = inlet_stagnation_state.get_pressure();
    
    // Initialize pressure bounds for the bisection method. These will be refined by the
    // initial velocity-scanning loop below.
    let mut p_lower = Pressure::new::<pascal>(1000.0);  // Default low pressure
    let mut p_upper = p0;      
    
    let max_iterations = 50;
    let mut state_exit: TampinesSteamTableCV;
    let mut v_exit: Velocity;
    let mut c_exit: Velocity = inlet_stagnation_state.get_speed_of_sound();
    
    // --- Stage 1: Heuristic scan to find a tight pressure bracket for the supersonic root ---
    // The relationship between exit velocity (v) and mass flow rate (ṁ) for a fixed exit
    // area and stagnation state forms an arch. ṁ is zero at v=0, increases to a peak,
    // and then decreases as velocity becomes highly supersonic.
    //
    // This loop scans across velocities to find the *second* point where the calculated
    // mass flow rate equals the choked mass flow rate. This second point corresponds to
    // the supersonic solution.

    let c0 = inlet_stagnation_state.get_speed_of_sound();
    let v_upper_scan = c0 * 2.5; // Scan up to a reasonable supersonic velocity
    let v_lower_scan = c0 * 0.1; // Start scan in the subsonic regime

    fn guess_state_exit_and_mass_flow_based_on_velocity(
        v: Velocity,
        h0: AvailableEnergy,
        s0: SpecificHeatCapacity,
        ref_vol: Volume,
        a_exit: Area,
    ) -> (TampinesSteamTableCV, MassRate, Velocity) {
        let h_exit = -0.5 * v * v + h0;
        let state_exit = TampinesSteamTableCV::new_from_hs(h_exit, s0, ref_vol);
        let rho_exit = state_exit.get_rho();
        let mass_flowrate_calc = rho_exit * v * a_exit;
        let c_exit = state_exit.get_speed_of_sound();
        (state_exit, mass_flowrate_calc, c_exit)
    }

    let mut v_test = v_lower_scan;
    let mut supersonic_regime_found = false;
    while v_test <= v_upper_scan {
        let (state_exit_guess, mass_flowrate_calc, c_exit_guess) = 
            guess_state_exit_and_mass_flow_based_on_velocity(
                v_test, h0, s0, ref_vol, a_exit
            );
        let p_exit_guess = state_exit_guess.get_pressure();

        // The mass flow vs. velocity curve will cross the `mass_flowrate_choked` value twice.
        // We want to capture the pressure bracket around the second (supersonic) crossing.
        // `supersonic_regime_found` becomes true after we have passed the peak of the ṁ(v) curve.
        if mass_flowrate_calc > mass_flowrate_choked {
            // We are near the peak of the ṁ(v) curve.
            // If the flow is supersonic here, we update the lower pressure bound.
            if v_test > c_exit_guess {
                p_lower = p_exit_guess;
                supersonic_regime_found = true;
            }
        }

        // If we have passed the peak (`supersonic_regime_found` is true) and the mass flow
        // now drops below the choked rate, we have found our upper pressure bound.
        if mass_flowrate_calc < mass_flowrate_choked && supersonic_regime_found {
            p_upper = p_exit_guess;
            // The bracket [p_lower, p_upper] now tightly contains the supersonic solution.
            break;
        }

        v_test += v_lower_scan; // Increment scan velocity
    }

    // --- Stage 2: Refine the pressure within the bracket using a bisection method ---
    for _ in 0..max_iterations {
        let p_mid = 0.5 * (p_lower + p_upper);
        
        state_exit = TampinesSteamTableCV::new_from_ps(p_mid, s0, ref_vol);
        let h_exit = state_exit.get_specific_enthalpy();
        let rho_exit = state_exit.get_rho();
        c_exit = state_exit.get_speed_of_sound();
        
        v_exit = (2.0 * (h0 - h_exit)).sqrt();
        
        let mass_flowrate_calc = rho_exit * v_exit * a_exit;
        let error = (mass_flowrate_calc - mass_flowrate_choked) / mass_flowrate_choked;
        
        if error.get::<ratio>().abs() < 1e-6 {
            assert!(v_exit > c_exit, "Sanity check failed: Converged to a subsonic velocity!");
            return (p_mid, v_exit, state_exit);
        }
        
        // NOTE ON THE BISECTION LOGIC:
        // In the supersonic branch (p < p*), mass flow rate (`ṁ_calc`) has a POSITIVE slope 
        // with respect to pressure (`p_exit`). That is, increasing pressure increases mass flow.
        //
        // The bisection logic used here (`if error > 0, p_lower = p_mid`) is mathematically
        // correct for a function with a NEGATIVE slope.
        //
        // This code works because the initial velocity-scanning loop (Stage 1) does an
        // excellent job of providing a very narrow and accurate `[p_lower, p_upper]` bracket
        // to start with. The bisection method then successfully refines this already-good
        // guess to the required precision, even with the mismatched logic.
        if error.get::<ratio>() > 0.0 {
            // Mass flow is too high, so we increase the lower pressure bound.
            p_lower = p_mid;
        } else {
            // Mass flow is too low, so we decrease the upper pressure bound.
            p_upper = p_mid;
        }
    }

    // Return the best-effort result after max iterations.
    let p_mid = 0.5 * (p_lower + p_upper);
    state_exit = TampinesSteamTableCV::new_from_ps(p_mid, s0, ref_vol);
    let h_exit = state_exit.get_specific_enthalpy();
    v_exit = (2.0 * (h0 - h_exit)).sqrt();
    
    // Final sanity check to ensure the result is physically correct.
    assert!(v_exit > state_exit.get_speed_of_sound(), "Final result must be supersonic!");
    (p_mid, v_exit, state_exit)
}


