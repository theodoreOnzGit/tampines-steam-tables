use uom::si::{f64::*, ratio::ratio, volume::cubic_meter};

use crate::prelude::{TampinesSteamTableCV, functional_programming::{hs_flash_eqm::{p_hs_eqm, w_hs_eqm}, ps_flash_eqm::v_ps_eqm}};


/// This is an algorithm to obtain outlet thermodynamic state 
/// for a converging nozzle with subsonic flow
#[inline]
pub fn get_choked_flow_state_for_nozzle_subsonic(
    p1: Pressure,
    h1: AvailableEnergy,
    a1: Area,
    a2: Area,
    v1: Velocity,
) -> (Pressure, AvailableEnergy, MassRate) {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);



    let stagnation_enthalpy: AvailableEnergy = 
        h1 + 0.5 * v1 * v1;

    let s1 = state_1.get_specific_entropy();

    // now, i'll have to get a solver for choked flow 
    // to obtain the outlet state, 
    // because speed of sound will actually depend on outlet state.
    // 
    // So I will need to iterate

    let mut c_guess: Velocity = state_1.get_speed_of_sound();
    // this is energy and mass balance
    let mut h2_guess: AvailableEnergy = 
        stagnation_enthalpy - 0.5 *c_guess * c_guess;

    let tolerance: f64 = 1e-5;
    let max_iter = 100;
    let mut n_iter = 0;

    // this will iterate until energy balance is okay, and speed of 
    // sound convergees
    loop {
        // obtain next guess of c_guess 
        let c_old = c_guess;
        c_guess = w_hs_eqm(h2_guess, s1);
        h2_guess = stagnation_enthalpy - 0.5 * c_guess * c_guess;

        // find residual of c 

        let c_residual: Ratio = (c_old- c_guess)/c_old;
        // if residual is less than some tolerance, break out
        if c_residual.abs() <= Ratio::new::<ratio>(tolerance) {
            break;
        }

        n_iter += 1;

        if n_iter >= max_iter {

            dbg!("maximum iterations reached for choked flow solver");
        }
        
    }

    // after this, it is time to do force balance
    // force balance is 
    //
    // P1A1 + v1 * mass_flowrate  = P2A2 + c*mass_flowrate
    let rho1: MassDensity = state_1.get_rho();
    let h2 = h2_guess;
    let s2 = s1;
    let c = c_guess;
    let p2 = p_hs_eqm(h2, s2);

    let p1a1: Force = p1 * a1;
    let p2a2: Force = p2 * a2;

    // for nozzle, 
    // mass flowrate =  (p1a1 - p2a2) / (v1 + c)

    let mass_flowrate = (p1a1 - p2a2)/(v1 + c);
    
    return (p2, h2, mass_flowrate);

}
/// this is for the special case of isentropy, 
/// but most of the time, 
/// it won't work for isentropy for choked flow scenario
/// as in mass and energy balance solves at the same time
#[inline]
fn force_balance_isentropic_nozzle(
    p1: Pressure,
    p2: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    ) -> Force {
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

    let rho1: MassDensity = state_1.get_specific_volume().recip();
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();
    // note: this is isentropic, hence s2 = s1;
    let s2 = s1;
    let v1: Velocity = mass_flowrate/rho1/a1;

    let p1a1: Force = p1*a1;
    // left hand side of momentum balance
    //
    // P1 A1 + dot{m}^2/{rho1 a1}
    let lhs: Force = p1a1 + mass_flowrate * v1;

    // if we have p2, we can get the second thermodynamic state
    let p2a2 = p2 * a2;
    let rho2 = v_ps_eqm(p2, s2).recip();
    let rhs: Force = p2a2 + mass_flowrate * mass_flowrate/rho2/a2;

    return lhs-rhs;

}
