use uom::si::{f64::*, volume::cubic_meter};

use crate::prelude::TampinesSteamTableCV;


/// This is an algorithm to obtain outlet thermodynamic state 
/// for a converging nozzle with subsonic flow
#[inline]
pub fn get_choked_flow_state_for_nozzle_subsonic(
    p1: Pressure,
    h1: AvailableEnergy,
    a1: Area,
    mass_flowrate: MassRate,
){

    let ref_vol = Volume::new::<cubic_meter>(1.0):
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let rho1: MassDensity = state_1.get_rho();

    let v1: Velocity = mass_flowrate/a1/rho1;


    let stagnation_enthalpy: AvailableEnergy = 
        h1 + 0.5 * v1 * v1;

    // now, i'll have to get a solver for choked flow 
    // to obtain the outlet state, 
    // because speed of sound will actually depend on outlet state.
    //

    let mut c_guess: Velocity = state_1.get_speed_of_sound();

    // for enthalpy, we 

    todo!();
    

}
