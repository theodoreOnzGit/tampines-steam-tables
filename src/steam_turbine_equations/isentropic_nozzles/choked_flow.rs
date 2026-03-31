use uom::si::{f64::*, ratio::ratio, volume::cubic_meter};

use crate::prelude::TampinesSteamTableCV;


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
pub fn get_choked_flow_state_for_nozzle_subsonic_to_sonic(
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
    use uom::si::area::square_centimeter;
    use uom::si::available_energy::kilojoule_per_kilogram;
    use uom::si::heat_capacity::kilojoule_per_kelvin;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::{kilopascal, megapascal};
    use uom::si::f64::*;
    use uom::si::ratio::ratio;
    use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::velocity::meter_per_second;
    use uom::si::volume::cubic_meter;

    use crate::prelude::TampinesSteamTableCV;
    use crate::steam_turbine_equations::choked_flow::get_choked_flow_state_for_nozzle_subsonic_to_sonic;


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
        let ref_vol = Volume::new::<cubic_meter>(1.0);

        let state_0 = TampinesSteamTableCV::new_from_tp_quality_1(
            t1, p0, ref_vol
        );

        let h0 = state_0.get_specific_enthalpy();

        // let's get entropy
        let s0 = state_0.get_specific_entropy();

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
            max_relative=5e-3,
        );

        // now for throat conditions 

        let state_throat = 
            TampinesSteamTableCV::new_from_ps(
                p_throat, s0, ref_vol
            );

        // let's find the velocity here 
        // this is, of course, the speed of sound
        let c = state_throat.get_speed_of_sound();

        // now, based on enthalpy balance, we should get 
        // a h_throat 

        let h_throat = state_throat.get_specific_enthalpy();

        approx::assert_relative_eq!(
            h_throat.get::<kilojoule_per_kilogram>(),
            3076.8,
            max_relative=1e-3,
        );

        let v_throat: Velocity = (2.0 *(h0-h_throat)).sqrt();

        // if we can see here, v_throat and c is about the same
        // Cengel's answer is 585.8 m/s
        // which is about there
        approx::assert_relative_eq!(
            c.get::<meter_per_second>(),
            584.74,
            max_relative=1e-3,
        );

        approx::assert_relative_eq!(
            v_throat.get::<meter_per_second>(),
            584.74,
            max_relative=1e-3,
        );

        // now, we find throat area is 10.33 cm^2 from Cengel's 
        // that's the answer, but we will work backwards using our function

        // velocity at the inlet is based on the mass flowrate 

        let a2 = Area::new::<square_centimeter>(10.33);
        // now, this a1 is arbitrary.., we only need know that v1 is negligble
        let a1 = 100.0 * a2;
        let rho_1 = state_0.get_rho();
        let v1: Velocity = mass_flowrate/rho_1/a1;
        
        approx::assert_relative_eq!(
            v1.get::<meter_per_second>(),
            3.659,
            max_relative=1e-3,
        );

        // let's see whether given these areas, we get the correct
        // critical pressure predictions 

        let (p_crit, h_crit, mass_flowrate) = 
            get_choked_flow_state_for_nozzle_subsonic_to_sonic(
                p1, h0, a1, a2, v1
            );

        approx::assert_relative_eq!(
            p_crit.get::<megapascal>(),
            1.092,
            max_relative=1e-3,
        );
        approx::assert_relative_eq!(
            h_crit.get::<kilojoule_per_kilogram>(),
            3076.8,
            max_relative=1e-3,
        );
        approx::assert_relative_eq!(
            mass_flowrate.get::<kilogram_per_second>(),
            2.5,
            max_relative=1e-3,
        );

        dbg!(&(
                critical_pressure_ratio,
                p_throat.get::<megapascal>(),
                v_throat.get::<meter_per_second>(),
                c.get::<meter_per_second>(),
        ));

        // now we want to deal with exit properties,
        // we know the pressure and 
        // Cengel writes the entahlpy as 2816.1 kJ/kg given a 0.93 
        // efficiency (we are not testing for that here)

        let h_exit = AvailableEnergy::new::<kilojoule_per_kilogram>(2816.1);

        let state_exit = 
            TampinesSteamTableCV::new_from_ph(
                p_exit, h_exit, ref_vol
            );

        let c_exit: Velocity = state_exit.get_speed_of_sound();


        // from Cengel's the exit speed of sound is 515.4 m/s
        approx::assert_relative_eq!(
            c_exit.get::<meter_per_second>(),
            515.4,
            max_relative=1e-3,
        );

        // from Cengel's the exit entropy is 7.2019 Kj/kg k

        let s_exit = state_exit.get_specific_entropy();

        approx::assert_relative_eq!(
            s_exit.get::<kilojoule_per_kilogram_kelvin>(),
            7.2019,
            max_relative=1e-3,
        );

        // for exit velocity, it is convenient to use a stagnation enthalpy 
        // as the reference, since it is the same through the 
        // whole nozzle 

        let v_exit = (2.0*(h0 - h_exit)).sqrt();
        approx::assert_relative_eq!(
            v_exit.get::<meter_per_second>(),
            929.8,
            max_relative=1e-3,
        );

        // exit mach number is 1.804 

        let mach_number_exit = state_exit.get_mach_number(v_exit);
        approx::assert_relative_eq!(
            mach_number_exit.get::<ratio>(),
            1.804,
            max_relative=1e-3,
        );

        dbg!(&(
                p_exit.get::<megapascal>(),
                v_throat.get::<meter_per_second>(),
                c_exit.get::<meter_per_second>(),
                v_exit.get::<meter_per_second>(),
                s_exit.get::<kilojoule_per_kilogram_kelvin>(),
                mach_number_exit.get::<ratio>(),
        ));

    }

}


