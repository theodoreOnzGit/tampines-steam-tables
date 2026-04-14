use uom::ConstZero;
use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::pressure::bar;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::{TampinesSteamTableCV};

// so here is the problem statement. 
//
//
// I have a converging diverging nozzle 
//
// there is a SET 
// inlet area, 
// throat area,
// outlet area,
//
// the conditions at inlet are known (this could be mass flowrate and 
// inlet velocity and such)
//
// The pressure at the turbine outlet will be known 
// (condenser pressure)
//
// It is very likely either 
//
// 1. saturated steam
// 2. superheated steam
//
// Water can in theory flow through the tubrine but that is highly 
// unlikely in steam turbine
//
// The turbine will have inlet pressure, and outlet pressure at condenser
//
// There will be a series of nozzles and turbine blades
// before then arriving at the outlet
// 
// The steam turbine can be a stator rotor pair interlaced with 
// many small control volumes in between each stator and rotor
//
// 
// for such a case, then there will be an inlet control volume with 
// known state, 
//
// then a stator (nozzle) and rotor pair
//
// For such a case, this is a pressure driven flow 
//
// the pressure difference between inlet and outlet will drive the flow
// It will be treated like a flow between two pressure points.
// Mass flowrate will be such that the pressures will equalise.
//
// For now, an impulse turbine stage can look like this 
// take two adjacent control volumes ,
//
//
// obtain pressure difference between them 
// knowing outlet pressure, I would then take this to be the nozzle 
// outlet pressure, 
//
// 1. I would then solve for the nozzle flow using the pressure differences
// 2. for the turbine blade, there is torque exerted, and the velocity 
// will be slowed down on the way out. This is based on the velocity 
// coefficient of the impulse turbine.
// 3. the outlet stagnation enthalpy is then calculated to be the state 
// of the steam entering the next control volume
//
// when pressure differences are given, we solve for mass flowrate 
// and velocity individually.
//
// I suppose we could do finite differencing. 
// In doing so, we may run into numerical errors such as Courant number.
// May be needed to do a local timestepping of sorts.
// Do we use rhoPimpleFoam algorithm? I myself am not sure.
//
// But the equations would be sort of an exponential type yeah... 
//
// not sure if I want to go down this rabbit hole 
// but nevermind.
//
// Anyway for a nozzle, with pressure differences, and no specified 
// mass flowrate 
//
// 1. assume stagnation in the inlet. (neglect KE in the inlet), otherwise 
// add it in
// 2. calculate throat velocity (possibly sonic), does this work?
// 3. generally speaking, nozzles just obey mass and energy balance.
// Energy balance is:
// v_throat = (2.0 * (h0 - h_throat).sqrt()
//
// assuming inlet stagnation states are known (p0, h0), we can 
// calculate entropy. This will also give us the inlet density of 
// steam.
//
// Nozzle is isentropic. (s1 = s0 = s_throat)
//
// area of entrance and throat are known 
//
// rho1 * a1 * v1 = rho2 * a2 * v2
//
// from this, we can obtain:
//
// rho2 * v2/v1 = rho1 * a1 / (a2);
//
// rho2 * v2/v1 is constant from mass balance
//
// What we can do, is to guess a mass flowrate.
//
// After guessing mass flowrate, we obtain v1 first 
// now we have the product of rho2 * v2.
//
// We have to guess h_throat iteratively. we shall have nested loop (oops)
//
// Otherwise, guess h_throat iteratively, the limits will be between 
// v_throat is supersonic and v_throat is 0
//
// from then we can guess v_throat. 
//
// Better yet, just guess v_throat (v2) first.
//
// This makes it easy to guess h_throat. 
// When h_throat is known, s_throat is known, thermodynamic state is fully 
// known. 
//
// we then get a formula for v2/v1. 
//
// From this, we get v1. Which is the mass flowrate.
//
// Now, when v_throat is guessed, the pressure at the throat is also 
// guessed.
//
// We have to see whether the pressure at the throat is good for the 
// pressure bounds given.
//
// Once having a throat pressure, we can obtain an equation to get 
// the mass flowrate in the outside part of the nozzle given the 
// nozzle diverging part.
//
// When the mass flowrate of the nozzle and supersonic diffuser match,
// given the pressure difference set, 
// then we have a mass flowrate solution in the nozzle.
//
//
// After solving the nozzle part, the mass flowrate is known. And the exit 
// velocity of the c-d nozzle is known. 
//
// Then we can just perform velocity reduction in the spinning part of 
// the impulse turbine. 
//
// And there we have it, we have mass flowrate and enthalpy at the end 
// of the turbine stage.
//
//
//
//


// given a sonic flow, 
//
// note, shocks may occur here 
// 
// given a pressure at the outlet, p2,
// and throat state, guess the state of flow going out
// mass flowrate is based on choked flow
//
// stagnation properties should also be supplied to facilitate calculation
//
// note that this is no longer isentropic
#[inline]
pub fn guess_state_for_diverge_nozzle_from_throat(
    p0: Pressure,
    h0: AvailableEnergy,
    p2: Pressure,
    a_throat: Area,
    a_out: Area,
    mass_rate_throat: MassRate,
    state_throat: TampinesSteamTableCV,
) -> TampinesSteamTableCV {
    // let's have a reference mass flux first 

    let mass_velocity_ref: MassFlux = mass_rate_throat/a_throat;

    // first I need throat pressure,

    let p_throat = state_throat.get_pressure();

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


    let mut h_outlet: AvailableEnergy;
    if (mass_velocity_error).abs() < 0.0001 {
        h_outlet = h2_ideal;
        return TampinesSteamTableCV::new_from_ph(
            p2, h_outlet, ref_vol
        );
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
            return TampinesSteamTableCV::new_from_ph(
                p2, h_outlet, ref_vol
            );
        }

    //    // Adjust bounds
    //    if mach_value < 1.0 {
    //        p_high = p_mid; // Need lower pressure (more expansion)
    //    } else {
    //        p_low = p_mid;  // Need higher pressure (less expansion)
    //    }

    //    // Check convergence
    //    if (p_high - p_low) < tolerance {
    //        return (p_low + p_high) / 2.0;
    //    }
    }

    // Return midpoint if not converged
    let h_mid_bound = 0.5 * (upper_bound + lower_bound);



    todo!()
}
// 
//
//
#[inline]
pub fn guess_massrate_and_state_for_converge_nozzle_from_stagnation(
    // these are stagnation pressure and enthalpy
    p0: Pressure,
    h0: AvailableEnergy,
    v_throat: Velocity,
    a_throat: Area,
    ) -> (MassRate, TampinesSteamTableCV) {
    // we have inlet stagnation enthalpy
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    // now we have inlet stagnation state
    let state_0 = TampinesSteamTableCV::new_from_ph(p0, h0, ref_vol);


    // calculate entropy, because process should be isentropic,
    // entropy is constant
    let s0 = state_0.get_specific_entropy();


    // calculate critical pressure
    let critical_pressure_ratio: Ratio = 
        state_0.get_critical_pressure_ratio();

    // this is critical pressure for mach 1
    let p_critical = critical_pressure_ratio * p0;
    let s_throat = s0;

    let state_choked = TampinesSteamTableCV::new_from_ps(
        p_critical, s_throat, ref_vol
    );

    // this is
    let mut state_out = state_choked;

    // we get speed of sound here 

    let c = state_choked.get_speed_of_sound();

    // this is the velocity for mass flowrate and energy balance 
    // calculation 
    let mut v = v_throat;

    // we know throat velocity cannot exceed mach 1
    if v >= c {
        // we are in a choked state, this is sonic flow
        v = c;
    } else {
        // we are not in a choked state, this is subsonic
        // in that case,
        let h_throat = h0 - 0.5 * v * v;
        state_out = TampinesSteamTableCV::new_from_hs(
            h_throat, s_throat, ref_vol
        );
        
    };
    let rho_out = state_out.get_rho();


    let mass_flowrate = rho_out * v * a_throat;




    // now, based on the above algorithm, we have to guess a mass 
    // flowrate given a pressure across the c-d nozzle 
    // but to do so, we have to iteratively guess velocity
    // until the mass flowrate in both converging and diverging parts 
    // of the nozzle tally
    
    return (mass_flowrate, state_out);

}


/// This gives the mass flowrate for choked flow given the area of a throat
/// and the inlet stagnation properties (negligble KE)
#[inline]
pub fn get_choked_flow_massrate_and_state_from_stagnation_properties_and_area(
    p0: Pressure,
    h0: AvailableEnergy,
    a_throat: Area,
) -> (MassRate, TampinesSteamTableCV) {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    // now we have inlet stagnation state
    let state_0 = TampinesSteamTableCV::new_from_ph(p0, h0, ref_vol);


    // calculate entropy, because process should be isentropic,
    // entropy is constant
    let s0 = state_0.get_specific_entropy();


    // calculate critical pressure
    let critical_pressure_ratio: Ratio = 
        state_0.get_critical_pressure_ratio();

    // this is critical pressure for mach 1
    let p_critical = critical_pressure_ratio * p0;
    let s_throat = s0;

    let state_choked = TampinesSteamTableCV::new_from_ps(
        p_critical, s_throat, ref_vol
    );

    let c = state_choked.get_speed_of_sound();

    let rho_choked = state_choked.get_rho();

    let massrate: MassRate = c * rho_choked * a_throat;


    return (massrate, state_choked);
}
