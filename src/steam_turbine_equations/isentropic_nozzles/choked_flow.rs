use uom::si::{f64::*, ratio::ratio, volume::cubic_meter};

use crate::prelude::{TampinesSteamTableCV, functional_programming::ps_flash_eqm::v_ps_eqm};


/// This is an algorithm to obtain outlet thermodynamic state 
/// for a converging nozzle with subsonic flow
/// Given inlet conditions (p1, h1, v1) and geometry (a1, a2), calculates
/// the throat conditions assuming choked flow (M = 1 at exit).
///
/// # Arguments
/// * `p1` - Inlet pressure
/// * `h1` - Inlet specific enthalpy
/// * `a1` - Inlet area
/// * `a2` - Throat (exit) area
/// * `v1` - Inlet velocity
///
/// # Returns
/// Tuple of (p2, h2, mass_flowrate) where:
/// * `p2` - Throat pressure
/// * `h2` - Throat specific enthalpy
/// * `mass_flowrate` - Choked mass flow rate (determined by throat conditions)
///
/// # Warnings
/// - Warns if mass balance error > 5% (inlet vs throat)
/// - Warns if momentum balance error > 5%
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
    let rho_1 = state_1.get_rho();


    let s1 = state_1.get_specific_entropy();

    let h0 = h1 + 0.5 * v1 * v1;
    // stagnation state (initial)
    let state_0 = TampinesSteamTableCV::new_from_hs(h0, s1, ref_vol);
    let p0 = state_0.get_pressure();


    // now, i'll have to get a solver for choked flow 

    // let's use the critical pressure 

    let critical_pressure_ratio: Ratio = 
        state_0.get_critical_pressure_ratio();

    // this is critical pressure for mach 1
    let p2 = critical_pressure_ratio * p0;
    // let's get speed of sound here 
    let s2 = s1;
    let state_2 = TampinesSteamTableCV::new_from_ps(p2, s2, ref_vol);
    let c = state_2.get_speed_of_sound();
    let h2 = state_2.get_specific_enthalpy();
    let rho_2 = state_2.get_rho();

    // now let's get mass balance first 

    let mass_flowrate = a1 * v1 * rho_1;
    let choked_mass_flowrate = a2 * c * rho_2;

    // assert that it isn't too different 
    //

    let mass_flowrate_error: Ratio = 
        (choked_mass_flowrate - mass_flowrate)/mass_flowrate;

    if mass_flowrate_error.abs() >= Ratio::new::<ratio>(0.05) {
        eprintln!("Warning: Mass balance error = {:.2}%. \
                   Inlet conditions may not be consistent with choked flow.",
                   mass_flowrate_error.get::<ratio>() * 100.0);
    }
    

    // after this, it is time to do force balance
    // force balance is 
    //
    // P1A1 + v1 * mass_flowrate  = P2A2 + c*mass_flowrate

    let p1a1: Force = p1 * a1;
    let p2a2: Force = p2 * a2;
    
    let lhs: Force = p1a1 + v1 * choked_mass_flowrate;
    let rhs: Force = p2a2 + c * choked_mass_flowrate;

    // for nozzle, 
    // mass flowrate =  (p1a1 - p2a2) / (-v1 + c)

    let momentum_balance_error: Ratio = (lhs - rhs)/lhs;

    if momentum_balance_error.abs() >= Ratio::new::<ratio>(0.05) {
        eprintln!("Warning: Momentum balance error = {:.2}%. \
                   Check if flow is truly choked.",
                   momentum_balance_error.get::<ratio>() * 100.0);
    }

    
    return (p2, h2, choked_mass_flowrate);

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
