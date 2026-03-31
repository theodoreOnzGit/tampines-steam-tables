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

#[cfg(test)]
mod choked_flow_examples{
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::{kilopascal, megapascal};
    use uom::si::f64::*;
    use uom::si::ratio::ratio;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::volume::cubic_meter;

    use crate::prelude::TampinesSteamTableCV;


    /// this is from example 17-16 in Cengel's thermodynamics 8th edition
    ///
    /// Steam enters a converging-diverging nozzle at 2 MPa and 400C 
    /// with negligble velocity and flowrate is 2.5 kg/s 
    /// it exits at pressure of 300 kPa 
    ///
    /// from nozzle entrance to throat, it is isentropic 
    /// after throat, it is 93 % efficient 
    ///
    /// what is 
    /// 1. the throat and exit area
    /// 2. mach number at throat and nozzle exit
    #[test]
    pub fn steam_flow_cd_nozzle(){

        let p1 = Pressure::new::<megapascal>(2.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(400.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(2.5);
        let p_exit = Pressure::new::<kilopascal>(300.0);

        // stagnation pressure = p1 as velocity is negligble 

        let p0 = p1;
        // we assume steam is superheated, so quality is 1 
        let x0 = 1.0;
        let ref_vol = Volume::new::<cubic_meter>(1.0);

        let state_0 = TampinesSteamTableCV::new_from_tp_quality_1(
            t1, p1, ref_vol
        );

        let critical_pressure_ratio = 
            state_0.get_critical_pressure_ratio();

        // In Cengel's textbook, critical pressure ratio is approximated 
        // as 0.546 

        approx::assert_relative_eq!(
            critical_pressure_ratio.get::<ratio>(),
            0.546,
            max_relative=1e-3,
        );

        // critical pressure ratio checks out!
        // now, we expect throat pressure to be 1.09 MPa
        let p_throat = p1 * critical_pressure_ratio;
        approx::assert_relative_eq!(
            p_throat.get::<megapascal>(),
            1.09,
            max_relative=1e-3,
        );



        dbg!(&(
                critical_pressure_ratio,
                p_throat,
        ));
        todo!();




    }

}


