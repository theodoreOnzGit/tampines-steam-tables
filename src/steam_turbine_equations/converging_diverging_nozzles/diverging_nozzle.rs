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
    // let's have a reference mass flux first 

    let mass_velocity_ref: MassFlux = mass_rate_throat/a_throat;

    // first I need throat pressure,

    // then we need to compare this to outlet pressure 
    //
    // if p_throat > p2, then good, we should continue having flow
    //
    // otherwise, we should expect deceleration 
    //

    // now at some threshold pressure, we should get smooth isentropic 
    // acceleration to supersonic flow
    //
    // This can be done using a (p,s) or (h,s) 
    // flash, and ensuring the mass balance 
    // holds

    let s2_ideal = state_throat.get_specific_entropy();
    // so we assume isentropy first

    let ref_vol = Volume::new::<cubic_meter>(1.0);


    let isentropic_outlet_state = 
        TampinesSteamTableCV::new_from_ps(
            p2, s2_ideal, ref_vol
        );

    let h2_ideal = isentropic_outlet_state.get_specific_enthalpy();


    // perhaps in general, a (p,h) algorithm may work... so long as the 
    // mass flowrate is satisfied
    //
    // this would abstract away any irreversibility

    fn guess_mass_velocity_given_ph_flash(
        h0: AvailableEnergy,
        p2: Pressure,
        h2_guess: AvailableEnergy) -> MassFlux {

        // note that this uses energy balance equations
        let v2: Velocity = (2.0 * (h0-h2_guess)).sqrt();
        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let outlet_state = 
            TampinesSteamTableCV::new_from_ph(
                p2, h2_guess, ref_vol
            );
        let rho2 = outlet_state.get_rho();
        let mass_velocity: MassFlux = v2 * rho2;

        return mass_velocity;
    }

    // the first case is the ideal case, if mass velocity is more 
    // than 
    let mass_velocity_ideal = 
        guess_mass_velocity_given_ph_flash(
            h0, p2, h2_ideal
        );
    let mass_velocity_error: f64 = 
        (
            (mass_velocity_ideal - mass_velocity_ref)/mass_velocity_ref
        ).get::<ratio>();


    let h_outlet: AvailableEnergy;
    if (mass_velocity_error).abs() < 0.0001 {
        h_outlet = h2_ideal;
        let v2: Velocity = (2.0 * (h0- h_outlet)).sqrt();
        return (v2,TampinesSteamTableCV::new_from_ph(
                p2, h_outlet, ref_vol
        ));
    }

    // what would be our upper and lower bounds for h2?
    // if we are going to experience shocks and subsonic flow, 
    // h2 > h_throat, but it will be lower than stagnation enthalpy 
    // so h0 should be the upper bound
    // The lower bound should be h2_isentropic, since this is the 
    // lowest possible given an isentropic state (higher h2 is irreversible)
    //
    let mut upper_bound = h0;
    let mut lower_bound = h2_ideal;




    let max_iterations = 50;
    // in these cases, we need to expect perhaps an accelerating 
    // section.
    //
    // Good videos to watch... 
    // Nozzle efficiency:
    // https://www.youtube.com/watch?v=qWq27t6sN6A
    //
    //
    // Here is a seires on supersonic and transonic cd nozzles
    // https://www.youtube.com/watch?v=GGrJXbkxRIs
    //

    // 1 kJ/kg tolerance
    let tolerance = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0); 

    for _ in 0..max_iterations {
        let h_mid_bound = 0.5 * (upper_bound + lower_bound);

        // Get mass velocity at this enthalpy 

        let h2_guess = h_mid_bound;

        let mass_velocity_guess = 
            guess_mass_velocity_given_ph_flash(
                h0, p2, h2_guess
            );

        let mass_velocity_error: f64 = 
            (
                (mass_velocity_guess - mass_velocity_ref)/mass_velocity_ref
            ).get::<ratio>();



        if (mass_velocity_error).abs() < 0.0001 {
            h_outlet = h2_guess;
            let v2: Velocity = (2.0 * (h0- h_outlet)).sqrt();
            return (v2,TampinesSteamTableCV::new_from_ph(
                p2, h_outlet, ref_vol
            ));
        }
        // if mass velocity error > 0 , we are guessing too high a mass 
        // flowrate, we may want to decrease enthalpy

        // Adjust bounds
        if mass_velocity_error > 0.0 {
            // mass velocity too high, need higher h2 (lower velocity)
            lower_bound = h_mid_bound;  
        } else {
            // mass velocity too low, need lower h2 (higher velocity)
            upper_bound = h_mid_bound; 
        }

        // Check convergence
        if (upper_bound - lower_bound) < tolerance {
            h_outlet = 0.5 * (upper_bound + lower_bound);
            let v2: Velocity = (2.0 * (h0- h_outlet)).sqrt();
            return (v2,TampinesSteamTableCV::new_from_ph(
                p2, h_outlet, ref_vol
            ));
        }
    }

    // Return midpoint if not converged
    h_outlet = 0.5 * (upper_bound + lower_bound);
    let v2: Velocity = (2.0 * (h0- h_outlet)).sqrt();

    return (v2,TampinesSteamTableCV::new_from_ph(
        p2, h_outlet, ref_vol
    ));
}

