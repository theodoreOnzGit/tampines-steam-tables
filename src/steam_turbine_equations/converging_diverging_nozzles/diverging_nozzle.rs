use uom::ConstZero;
use uom::si::f64::*;
use uom::si::pressure::pascal;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::TampinesSteamTableCV;
use crate::steam_turbine_equations::joule_thomson::get_outlet_velocity_and_state_joule_thomson;

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
    a_exit: Area,
    mass_rate_throat: MassRate,
    state_throat: TampinesSteamTableCV,
) -> (Velocity, TampinesSteamTableCV) {
    
    // Calculate reference mass flowrate
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let inlet_stagnation_state = 
        TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);

    let p0 = inlet_stagnation_state.get_pressure();
    // Calculate perfectly expanded solution (supersonic branch)
    let (p_ideal_exp_supersonic, v_ideal_exp_supersonic, state_ideal_exp_supersonic) 
        = calculate_isentropic_exit_pressure_velocity_and_state_supersonic(
            inlet_stagnation_state, 
            a_exit, 
            mass_rate_throat,
        );

    // Calculate isentropic subsonic branch
    let (p_ideal_exp_subsonic, v_ideal_exp_subsonic, state_ideal_exp_subsonic) 
        = calculate_isentropic_exit_pressure_velocity_and_state_subsonic(
            inlet_stagnation_state, 
            a_exit, 
            mass_rate_throat,
        );
    let debug = false;

    if debug {
        dbg!(&(p_ideal_exp_subsonic,p_ideal_exp_supersonic));
        dbg!(&(v_ideal_exp_subsonic,v_ideal_exp_supersonic));
    }

    // Helper: calculate mass flowrate using outlet enthalpy (p,h) flash 
    // using velocit as input
    fn calculate_mass_rate_and_state_at_outlet_ph_velocity(
        h0: AvailableEnergy,
        p2: Pressure,
        v2: Velocity,
        a2: Area,
    ) -> (MassRate, TampinesSteamTableCV) {
        // Energy equation: v₂ = √(2(h₀ - h₂))
        //
        // we use: h2 = h0 - 0.5 * v2^2
        let h2: AvailableEnergy = h0 - 0.5 * v2 * v2;
        
        // Get density from (p,h) flash
        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_2 = TampinesSteamTableCV::new_from_ph(p2, h2, ref_vol);
        let rho2 = state_2.get_rho();
        
        // Mass flux: G = ρv
        // Mass rate: G*a2
        let mass_rate = rho2 * v2 * a2;
        
        (mass_rate, state_2)
    }

    let p_crit = state_throat.get_pressure();

    // now let's have a sanity check 
    // p0 > p_ideal_exp_subsonic > p_crit > p_ideal_exp_supersonic

    assert!(p0 > p_ideal_exp_subsonic);
    assert!(p_ideal_exp_subsonic > p_crit);
    assert!(p_crit > p_ideal_exp_supersonic);



    // so before anything, we have a few pressures to take note of 
    // 
    // Since we already assume choked flow
    // 
    // 1. stagnation pressure (the absolute upper bound)
    // 2. critical pressure (the pressure at the throat)
    // 3. subsonic ideal expansion pressure, that is if we happen to hv choked 
    // flow, and then it falls back down to subsonic speed isentropically,
    // then it is this pressure 
    //
    // 4. supersonic ideal expansion pressure, that is, if we happen to 
    // have choked flow, and then it accelerates to supersonic speed 
    // isentropically, it will leave the exit at this pressure
    //
    // 5. p2, the outlet pressure

    // first off, if pressure is higher than subsonic ideal expansion 
    // pressure, we should not even have choked flow in the first place 

    if p2 > p_ideal_exp_subsonic {
        eprintln!("outlet pressure is too high for choked flow to happen");
        eprintln!("outlet_pressure:");
        dbg!(&(p2));
        eprintln!("whereas isentropic subsonic expansion pressure is:");
        dbg!(&(p_ideal_exp_subsonic));
        panic!("");
    }


    // ========================================================================
    // Step 1: Try isentropic solution (no shocks)
    // ========================================================================

    // firstly, we want to check the nozzle boundary pressure,
    // it should be at least, the ideal expansion pressure
    // if p2 is lower than this ideal expansion pressure, oblique shocks 
    // will form outside
    // to save us some trouble, we will do this to within a tolerance 
    // of 10 Pa

    // or appropriate tolerance
    let pressure_tolerance = Pressure::new::<pascal>(10.0); 
    let pressure_diff_subsonic = (p2 - p_ideal_exp_subsonic).abs();
    let pressure_diff_supersonic = (p2 - p_ideal_exp_supersonic).abs();


    if pressure_diff_subsonic < pressure_tolerance {
        return (v_ideal_exp_supersonic, state_ideal_exp_supersonic);
    }

    if pressure_diff_supersonic < pressure_tolerance {
        return (v_ideal_exp_subsonic, state_ideal_exp_subsonic);
    }

    // ========================================================================
    // Step 2: Non-isentropic solution (shocks present) - Use Regula Falsi
    // ========================================================================
    //
    // in this case, p2 lies between the supersonic and subsonic 
    // expansion pressure branches, 
    // we should expect normal shocks
    //
    // So, the outlet thermodynamic state is fixed using p2
    // and we have a fixed (choked) mass flowrate
    // entropy won't be constant, so we cannot use that
    // 
    // we can vary velocity, 
    // calculate an enthalpy, 
    //
    // check the mass flowrate, and then guess the outlet state 
    //
    // for this, the velocity bounds should be between 
    // the supersonic ideal velocity and the subsonic velocity
    //
    // Shocks will occur, but the resulting velocities may or may not 
    // be supersonic.
    //
    // The bottom line is that mass flowrate must be conserved

    let max_iterations = 50;
    const TOLERANCE: f64 = 0.0001;  // 0.01% tolerance
    if p2 > p_ideal_exp_supersonic && p2 < p_ideal_exp_subsonic {

        // we are going to do a velocity scan algorithm again
        let mut v_upper_limit = v_ideal_exp_supersonic;
        let mut v_lower_limit = v_ideal_exp_subsonic;
        let v_increment = (v_upper_limit - v_lower_limit) * 0.1;

        // remember, we are supposed to vary v until the mass flowrate 
        // calculated reaches that of the throat

        let root_finder = |v_test: Velocity| -> (MassRate, TampinesSteamTableCV) {

            let (mass_flowrate_calc, outlet_state_with_shocks) = 
                calculate_mass_rate_and_state_at_outlet_ph_velocity(
                h0, p2, v_test, a_exit);

            let error = mass_flowrate_calc - mass_rate_throat;


            (error, outlet_state_with_shocks)
        };

        let mut v_test = v_lower_limit;
        let mut relative_error: f64 = 20_f64 * TOLERANCE;
        let (initial_error, outlet_state) = root_finder(v_lower_limit);

        // get the sign of the initial error

        let initial_error_positive: bool;

        if initial_error > MassRate::ZERO {
            initial_error_positive = true;
        } else {
            initial_error_positive = false;
        }


        if debug {
            dbg!(&(v_test));
            dbg!(&(initial_error,outlet_state));
        }
        

        // this is the initial velocity scan
        while v_test < v_upper_limit {
            let (test_error, outlet_state) = root_finder(v_test);

            // next is to check the sign of the error

            let test_error_positive: bool;

            if test_error > MassRate::ZERO {
                test_error_positive = true;
            } else {
                test_error_positive = false;
            }

            // if the signs are same, then continue 

            if initial_error_positive == test_error_positive {
                v_lower_limit = v_test;
                v_test += v_increment;

                if debug {
                    dbg!(&(v_test));
                    dbg!(&(test_error,outlet_state));
                }
                continue;
            };

            // if signs are not same

            if debug {
                dbg!(&(v_test));
                dbg!(&(test_error,outlet_state));
            }
            v_upper_limit = v_test; 
            break;

        }
        // now i can do bisection (or a secant method) 
        // between these two limits
        // since it's quite near the root
        //
        // or as AI suggested, I'm going to try Regula Falsi
        // near this region

        if debug{
            println!("Regula Falsi bounds found");
            dbg!(&(v_lower_limit,v_upper_limit));
        }

        let mut iteration = 0;

        let (mut error_lower_limit, _state_lower_limit) = 
            root_finder(v_lower_limit);
        let (mut error_upper_limit, _state_upper_limit) = 
            root_finder(v_upper_limit);

        if error_lower_limit.value * error_upper_limit.value >= 0.0 {
            panic!("bounds are same sign!");
        }

        // this is regula falsi
        while relative_error.abs() > TOLERANCE && iteration < max_iterations {

            // using secant formula
            v_test = 
                v_upper_limit - (v_upper_limit - v_lower_limit) * 
                error_upper_limit/(error_upper_limit - error_lower_limit);

            // check mass flowrate using velocity
            let (test_error, test_outlet_state) = root_finder(v_test);
            // update relative error
            relative_error = (test_error/mass_rate_throat).get::<ratio>();

            if relative_error.abs() < TOLERANCE {
                return (v_test, test_outlet_state);
            }

            // update bounds, keep root bracketed 

            if error_lower_limit.value * test_error.value < 0.0 {

                v_upper_limit = v_test;
                error_upper_limit = test_error;
            } else {

                v_lower_limit = v_test;
                error_lower_limit = test_error;
            }


            if debug {
                dbg!(&(v_lower_limit,v_upper_limit));
            }
            iteration += 1;

        }
        
        

    }

    // now, third part, is where the outlet pressure is lower than 
    // the perfectly expanded supersonic pressure.
    //
    // We expect shock waves to form outside the nozzle
    // This is beyond the scope of Cengel 
    //
    // This is in an underexpanded nozzle
    //
    // Now, for this to work, mass flowrate is the same. 
    // Moreover we know the nozzle pressure and eventual pressure 
    //
    // m_nozzle = m_out 
    //
    // p_nozzle known (p_ideal_exp_supersonic) --> This thermodynamic 
    // state is known
    // p_out known (p2)
    //
    // p_ideal_exp_supersonic > p2
    //
    // If area is constant, we can assume a joule thompson effect
    // I have this coded in another part of code

    if p2 < p_ideal_exp_supersonic {
        // joule thomson effect 

        let p1 = p_ideal_exp_supersonic;
        let h1 = state_ideal_exp_supersonic.get_specific_enthalpy();
        let (v2, state2) = 
            get_outlet_velocity_and_state_joule_thomson(
                p1, h1, p2, mass_rate_throat, a_exit
            );

        return (v2, state2);

    }
    
    // if none of these cases fit, panic

    panic!("Choked flow solver could not converge");

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
    let mut c_exit: Velocity;
    
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



#[inline]
pub fn calculate_isentropic_exit_pressure_velocity_and_state_subsonic(
    inlet_stagnation_state: TampinesSteamTableCV,
    a_exit: Area,
    mass_flowrate_choked: MassRate,
) -> (Pressure, Velocity, TampinesSteamTableCV) {
    
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    
    let h0: AvailableEnergy = inlet_stagnation_state.get_specific_enthalpy();
    let s0: SpecificHeatCapacity = inlet_stagnation_state.get_specific_entropy();
    let p0: Pressure = inlet_stagnation_state.get_pressure();
    
    // --- Stage 1: Set pressure bounds based on physics ---
    // The subsonic solution for a diverging nozzle must lie between the critical pressure
    // at the throat (p*) and the stagnation pressure (p0).
    let p_critical = inlet_stagnation_state.get_critical_pressure_ratio_pure_vapour() * p0;
    let mut p_lower = p_critical;
    let mut p_upper = p0;      
    
    let max_iterations = 50;
    let mut state_exit: TampinesSteamTableCV;
    let mut v_exit: Velocity;
    let mut c_exit: Velocity;

    // --- Stage 2: Refine the pressure using a bisection method ---
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
            assert!(v_exit < c_exit, "Sanity check failed: Converged to a supersonic velocity!");
            return (p_mid, v_exit, state_exit);
        }
        
        // Bisection logic for the SUBSONIC branch (negative slope):
        // In this regime, increasing pressure decreases mass flow rate.
        if error.get::<ratio>() > 0.0 {
            // Mass flow is too high. To reduce it, we must INCREASE the pressure.
            p_lower = p_mid;
        } else {
            // Mass flow is too low. To increase it, we must DECREASE the pressure.
            p_upper = p_mid;
        }
    }

    // Return the best-effort result after max iterations.
    let p_mid = 0.5 * (p_lower + p_upper);
    state_exit = TampinesSteamTableCV::new_from_ps(p_mid, s0, ref_vol);
    let h_exit = state_exit.get_specific_enthalpy();
    v_exit = (2.0 * (h0 - h_exit)).sqrt();
    
    // Final sanity check to ensure the result is physically correct.
    assert!(v_exit < state_exit.get_speed_of_sound(), "Final result must be subsonic!");
    (p_mid, v_exit, state_exit)
}
