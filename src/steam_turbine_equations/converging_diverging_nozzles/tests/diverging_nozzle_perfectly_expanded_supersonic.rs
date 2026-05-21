use uom::si::f64::*;
use uom::si::length::meter;
use uom::si::mass_rate::kilogram_per_second;
use uom::si::pressure::{bar, pascal};
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;
use uom::si::volume::cubic_meter;

use crate::prelude::{TampinesSteamTableCV, get_choked_flow_massrate_and_state_from_stagnation_properties_and_area};
use crate::steam_turbine_equations::diverging_nozzle::calculate_isentropic_exit_pressure_velocity_and_state_supersonic;
use crate::steam_turbine_equations::diverging_nozzle::calculate_isentropic_exit_pressure_velocity_and_state_subsonic;

/// this test checks the function for perfectly expanded dry steam
#[test]
fn diverging_nozzle_perfectly_expanded_supersonic_dry_steam(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let temperature = ThermodynamicTemperature::new::<degree_celsius>(400.0);
    let p1 = Pressure::new::<bar>(20.0);
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality_1(
        temperature, p1, ref_vol
    );

    let h1 = inlet_state.get_specific_enthalpy();
    let v1 = Velocity::new::<meter_per_second>(0.5);


    // now I'm going to use moore nozzle heights 
    // these were given by Google's AI, and Claude Sonnet
    // just need to check later
    let inlet_height = Length::new::<meter>(0.05635);
    let throat_height = Length::new::<meter>(0.05000);
    let exit_height = Length::new::<meter>(0.07200);
    let width = Length::new::<meter>(0.1);

    let a_throat = width * throat_height;
    let a2 = width * exit_height;
    let _a1 = width * inlet_height;

    // this is the stagnation state part
    let h0: AvailableEnergy = h1 + 0.5 * v1 * v1;
    // we get the entropy 
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();

    // now stagnation properties come from s1 = s0 
    let s0 = s1;
    let inlet_stagnation_state = TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);
    let p0 = inlet_stagnation_state.get_pressure();
    
    let (choked_mass_flowrate, choked_state) = 
        get_choked_flow_massrate_and_state_from_stagnation_properties_and_area(
            p0, h0, a_throat
        );
    let p_throat_critical = choked_state.get_pressure();

    // then now, we have the perfectly expanded part 


    let (p_exit_ideal, v_exit, state_exit) = 
        calculate_isentropic_exit_pressure_velocity_and_state_supersonic(
            inlet_stagnation_state, 
            a2, 
            choked_mass_flowrate
        );

    let exit_mass_flowrate: MassRate = 
        state_exit.get_rho()
        * v_exit 
        * a2;
    
    // now let's test the exit mass flowrate first

    approx::assert_relative_eq!(
        choked_mass_flowrate.get::<kilogram_per_second>(),
        exit_mass_flowrate.get::<kilogram_per_second>(),
        max_relative=1e-6
        );

    // next, we should assert that the exit state is isentropic

    approx::assert_relative_eq!(
        state_exit.get_specific_entropy().get::<joule_per_kilogram_kelvin>(),
        s0.get::<joule_per_kilogram_kelvin>(),
        max_relative=1e-6
        );

    // next, assert some sanity checks for pressure 
    //
    // the ideal exit pressure should be less than the throat pressure 
    //
    // and also we should have supersonic speed
    assert!(p_throat_critical > p_exit_ideal);
    let c_exit = state_exit.get_speed_of_sound();

    assert!(v_exit > c_exit);


    // for regression, let's see the exit pressure
    // should be about 379.2 kPa
    approx::assert_relative_eq!(
        p_exit_ideal.get::<pascal>(),
        379208_f64,
        max_relative=1e-6
        );

    // also want to check the mach number
    // is about 1.76 (quite reasonable)
    approx::assert_relative_eq!(
        state_exit.get_mach_number(v_exit).get::<ratio>(),
        1.76125,
        max_relative=1e-6
        );

}

/// this test checks the function for choked flow going back 
/// isentropically to subsonic flow
#[test]
fn diverging_nozzle_perfectly_expanded_subsonic_dry_steam(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let temperature = ThermodynamicTemperature::new::<degree_celsius>(400.0);
    let p1 = Pressure::new::<bar>(20.0);
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality_1(
        temperature, p1, ref_vol
    );

    let h1 = inlet_state.get_specific_enthalpy();
    let v1 = Velocity::new::<meter_per_second>(0.5);


    // now I'm going to use moore nozzle heights 
    // these were given by Google's AI, and Claude Sonnet
    // just need to check later
    let inlet_height = Length::new::<meter>(0.05635);
    let throat_height = Length::new::<meter>(0.05000);
    let exit_height = Length::new::<meter>(0.07200);
    let width = Length::new::<meter>(0.1);

    let a_throat = width * throat_height;
    let a2 = width * exit_height;
    let _a1 = width * inlet_height;

    // this is the stagnation state part
    let h0: AvailableEnergy = h1 + 0.5 * v1 * v1;
    // we get the entropy 
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();

    // now stagnation properties come from s1 = s0 
    let s0 = s1;
    let inlet_stagnation_state = TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);
    let p0 = inlet_stagnation_state.get_pressure();
    
    let (choked_mass_flowrate, choked_state) = 
        get_choked_flow_massrate_and_state_from_stagnation_properties_and_area(
            p0, h0, a_throat
        );
    let p_throat_critical = choked_state.get_pressure();

    // then now, we have the perfectly expanded part 


    let (p_exit_ideal, v_exit, state_exit) = 
        calculate_isentropic_exit_pressure_velocity_and_state_subsonic(
            inlet_stagnation_state, 
            a2, 
            choked_mass_flowrate
        );

    let exit_mass_flowrate: MassRate = 
        state_exit.get_rho()
        * v_exit 
        * a2;
    
    // now let's test the exit mass flowrate first

    approx::assert_relative_eq!(
        choked_mass_flowrate.get::<kilogram_per_second>(),
        exit_mass_flowrate.get::<kilogram_per_second>(),
        max_relative=1e-6
        );

    // next, we should assert that the exit state is isentropic

    approx::assert_relative_eq!(
        state_exit.get_specific_entropy().get::<joule_per_kilogram_kelvin>(),
        s0.get::<joule_per_kilogram_kelvin>(),
        max_relative=1e-6
        );

    // next, assert some sanity checks for pressure 
    //
    // the ideal exit pressure should be less than the throat pressure 
    //
    // and also we should have supersonic speed
    assert!(p_throat_critical < p_exit_ideal);
    let c_exit = state_exit.get_speed_of_sound();

    assert!(v_exit < c_exit);


    // for regression, let's see the exit pressure
    // should be about 1750 kPa
    approx::assert_relative_eq!(
        p_exit_ideal.get::<pascal>(),
        1750063_f64,
        max_relative=1e-6
        );

    // also want to check the mach number
    // is about 0.458 (quite reasonable)
    // below mach 1, subsonic is correct
    approx::assert_relative_eq!(
        state_exit.get_mach_number(v_exit).get::<ratio>(),
        0.4585318,
        max_relative=1e-6
        );

}

/// this test checks the function for perfectly expanded wet steam
#[test]
#[ignore="test not ready"]
fn diverging_nozzle_perfectly_expanded_supersonic_wet_steam(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let temperature = ThermodynamicTemperature::new::<degree_celsius>(400.0);
    let p1 = Pressure::new::<bar>(20.0);
    let quality = 0.2;
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality(
        temperature, p1, ref_vol,
        quality
    );

    let h1 = inlet_state.get_specific_enthalpy();
    let v1 = Velocity::new::<meter_per_second>(0.5);


    // now I'm going to use moore nozzle heights 
    // these were given by Google's AI, and Claude Sonnet
    // just need to check later
    let inlet_height = Length::new::<meter>(0.05635);
    let throat_height = Length::new::<meter>(0.05000);
    let exit_height = Length::new::<meter>(0.07200);
    let width = Length::new::<meter>(0.1);

    let a_throat = width * throat_height;
    let a2 = width * exit_height;
    let _a1 = width * inlet_height;

    // this is the stagnation state part
    let h0: AvailableEnergy = h1 + 0.5 * v1 * v1;
    // we get the entropy 
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();

    // now stagnation properties come from s1 = s0 
    let s0 = s1;
    let inlet_stagnation_state = TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);
    let p0 = inlet_stagnation_state.get_pressure();
    
    let (choked_mass_flowrate, choked_state) = 
        get_choked_flow_massrate_and_state_from_stagnation_properties_and_area(
            p0, h0, a_throat
        );
    let p_throat_critical = choked_state.get_pressure();

    // then now, we have the perfectly expanded part 


    let (p_exit_ideal, v_exit, state_exit) = 
        calculate_isentropic_exit_pressure_velocity_and_state_supersonic(
            inlet_stagnation_state, 
            a2, 
            choked_mass_flowrate
        );

    let exit_mass_flowrate: MassRate = 
        state_exit.get_rho()
        * v_exit 
        * a2;
    
    // now let's test the exit mass flowrate first

    approx::assert_relative_eq!(
        choked_mass_flowrate.get::<kilogram_per_second>(),
        exit_mass_flowrate.get::<kilogram_per_second>(),
        max_relative=1e-6
        );

    // next, we should assert that the exit state is isentropic

    approx::assert_relative_eq!(
        state_exit.get_specific_entropy().get::<joule_per_kilogram_kelvin>(),
        s0.get::<joule_per_kilogram_kelvin>(),
        max_relative=1e-6
        );

    // next, assert some sanity checks for pressure 
    //
    // the ideal exit pressure should be less than the throat pressure 
    //
    // and also we should have supersonic speed
    assert!(p_throat_critical > p_exit_ideal);
    let c_exit = state_exit.get_speed_of_sound();

    assert!(v_exit > c_exit);


    // for regression, let's see the exit pressure
    // should be about 379.2 kPa
    approx::assert_relative_eq!(
        p_exit_ideal.get::<pascal>(),
        379208_f64,
        max_relative=1e-6
        );

    // also want to check the mach number
    // is about 1.76 (quite reasonable)
    approx::assert_relative_eq!(
        state_exit.get_mach_number(v_exit).get::<ratio>(),
        1.76125,
        max_relative=1e-6
        );
    let inlet_state_quality_1 = TampinesSteamTableCV::new_from_tp_quality_1(
        temperature, p1, ref_vol
    );

    dbg!(&inlet_state_quality_1);
    dbg!(&inlet_state);

    // assert quality of inlet state
    approx::assert_relative_eq!(
        inlet_state.get_quality(),
        0.2,
        max_relative=1e-6
        );

    // check the quality of the exit state
    // (it should be 1 as the steam is fully vapour)
    approx::assert_relative_eq!(
        state_exit.get_quality(),
        1.0,
        max_relative=1e-6
        );
}
