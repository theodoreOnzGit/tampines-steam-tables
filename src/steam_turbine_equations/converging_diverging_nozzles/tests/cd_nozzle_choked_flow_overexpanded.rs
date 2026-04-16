use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::length::meter;
use uom::si::mass_rate::kilogram_per_second;
use uom::si::pressure::{bar, pascal};
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;
use uom::si::volume::cubic_meter;

use crate::prelude::TampinesSteamTableCV;
use crate::steam_turbine_equations::calculate_velocity_mass_flowrate_and_state_in_cd_nozzle;
// note: From google AI,
//
// International Test Series on Steam Nozzles
// is a good place to look for steam nozzle (validation)
// Moore Nozzles (specifically Nozzle B)
// IWSEP Nozzle
//
// The other thing:
// NASA CDV Nozzle Reference
//
// For moore nozzle B,
// throat area is 
//
// These are AI generated test cases
// Test Cases Created
//1. Subsonic Flow (No Choking)
//
//    Dry Steam: 10 bar, 300°C → 8 bar
//    Wet Steam: 10 bar, quality 0.95 → 8 bar
//    Checks: Isentropic flow, energy conservation, pressure match
//
//2. Choked Flow Back to Subsonic
//
//    Skipped (physically unlikely in CD nozzles - would need special geometry)
//
//3. Over-Expanded (Normal Shock Inside)
//
//    Dry Steam: 20 bar, 400°C → 12 bar (with shock)
//    Wet Steam: 15 bar, saturated → 10 bar
//    Checks: Entropy increases, non-isentropic, energy conserved
//
//4. Perfectly Expanded (Isentropic Throughout)
//
//    Dry Steam: 30 bar, 450°C → design pressure
//    Wet Steam: 20 bar, saturated → design pressure
//    Checks: Zero entropy change, matches ideal expansion
//
//5. Under-Expanded (Oblique Shocks Outside)
//
//    Dry Steam: 40 bar, 500°C → 60% of design pressure
//    Wet Steam: 25 bar, saturated → 50% of design pressure
//    Checks: Joule-Thomson throttling, entropy increases
//
#[test]
fn dry_steam_test(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let temperature = ThermodynamicTemperature::new::<degree_celsius>(400.0);
    let p1 = Pressure::new::<bar>(20.0);
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality_1(
        temperature, p1, ref_vol
    );
    let p2 = Pressure::new::<bar>(8.0);

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



    let (v, mass_rate, outlet_state) = 
        calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
            p1, 
            h1, 
            v1, 
            a_throat, 
            a2, 
            p2
        );


    // we are going to check for 
    // (1) isentropic flow 
    // (2) mass balance 
    // (3) energy balance 
    // (4) pressure match
    
    let s2 = outlet_state.get_specific_entropy();
    let s1 = inlet_state.get_specific_entropy();
    
    // we also consider stagnation state
    let h0_ref = h1 + 0.5 * v1 * v1;
    let state_0 = TampinesSteamTableCV::new_from_hs(h0_ref, s1, ref_vol);
    let s0 = state_0.get_specific_entropy();

    // For subsonic isentropic flow: s2 = s0 = s1
    println!("✓ Entropy Increased: s1 = {:.4} kJ/kg·K, s2 = {:.4} kJ/kg·K", 
        s1.get::<kilojoule_per_kilogram_kelvin>(),
        s2.get::<kilojoule_per_kilogram_kelvin>()
    );
    assert!(s2 > s1);
    assert!(s2 > s0);
    
    // ====================================================================
    // Test (2): Mass Balance - ṁ = ρ₂ v₂ A₂
    // ====================================================================
    let rho2 = outlet_state.get_rho();
    let mass_rate_calculated = rho2 * v * a2;

    approx::assert_relative_eq!(
        mass_rate.get::<kilogram_per_second>(),
        mass_rate_calculated.get::<kilogram_per_second>(),
        epsilon = 1e-6
    );

    println!("✓ Mass balance: ṁ = {:.6} kg/s", 
        mass_rate.get::<kilogram_per_second>()
    );
    // ====================================================================
    // Test (3): Energy Balance - h₀ = h₂ + v₂²/2
    // ====================================================================
    let h2 = outlet_state.get_specific_enthalpy();
    let h0_actual = h2 + 0.5 * v * v;

    approx::assert_relative_eq!(
        h0_ref.get::<kilojoule_per_kilogram>(),
        h0_actual.get::<kilojoule_per_kilogram>(),
        epsilon = 0.1  // 0.1 kJ/kg tolerance
    );

    println!("✓ Energy conserved: h₀ = {:.2} kJ/kg, h₂ + v₂²/2 = {:.2} kJ/kg", 
        h0_ref.get::<kilojoule_per_kilogram>(),
        h0_actual.get::<kilojoule_per_kilogram>()
    );
    // ====================================================================
    // Test (4): Pressure Match
    // ====================================================================
    let p2_result = outlet_state.get_pressure();

    approx::assert_relative_eq!(
        p2.get::<pascal>(),
        p2_result.get::<pascal>(),
        epsilon = 1e-3  // Very tight tolerance for pressure
    );

    println!("✓ Pressure match: p₂ = {:.2} bar", 
        p2_result.get::<bar>()
    );

}


#[test]
fn wet_steam_test(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let p1 = Pressure::new::<bar>(15.0);
    let x: f64 = 0.80;
    let inlet_state = TampinesSteamTableCV::new_from_sat_pressure_quality(
        p1, x, ref_vol
    );

    let p2 = Pressure::new::<bar>(10.0);

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



    let (v, mass_rate, outlet_state) = 
        calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
            p1, 
            h1, 
            v1, 
            a_throat, 
            a2, 
            p2
        );

    // we are going to check for 
    // (1) isentropic flow 
    // (2) mass balance 
    // (3) energy balance 
    // (4) pressure match
    
    let s2 = outlet_state.get_specific_entropy();
    let s1 = inlet_state.get_specific_entropy();
    
    // we also consider stagnation state
    let h0_ref = h1 + 0.5 * v1 * v1;
    let state_0 = TampinesSteamTableCV::new_from_hs(h0_ref, s1, ref_vol);
    let s0 = state_0.get_specific_entropy();

    // For subsonic isentropic flow: s2 = s0 = s1
    println!("s1 = {:.4} kJ/kg·K, s2 = {:.4} kJ/kg·K", 
        s1.get::<kilojoule_per_kilogram_kelvin>(),
        s2.get::<kilojoule_per_kilogram_kelvin>()
    );
    assert!(s2 > s1);
    assert!(s2 > s0);

    
    // ====================================================================
    // Test (2): Mass Balance - ṁ = ρ₂ v₂ A₂
    // ====================================================================
    let rho2 = outlet_state.get_rho();
    let mass_rate_calculated = rho2 * v * a2;

    approx::assert_relative_eq!(
        mass_rate.get::<kilogram_per_second>(),
        mass_rate_calculated.get::<kilogram_per_second>(),
        epsilon = 1e-6
    );

    println!("✓ Mass balance: ṁ = {:.6} kg/s", 
        mass_rate.get::<kilogram_per_second>()
    );
    // ====================================================================
    // Test (3): Energy Balance - h₀ = h₂ + v₂²/2
    // ====================================================================
    let h2 = outlet_state.get_specific_enthalpy();
    let h0_actual = h2 + 0.5 * v * v;

    approx::assert_relative_eq!(
        h0_ref.get::<kilojoule_per_kilogram>(),
        h0_actual.get::<kilojoule_per_kilogram>(),
        epsilon = 0.1  // 0.1 kJ/kg tolerance
    );

    println!("✓ Energy conserved: h₀ = {:.2} kJ/kg, h₂ + v₂²/2 = {:.2} kJ/kg", 
        h0_ref.get::<kilojoule_per_kilogram>(),
        h0_actual.get::<kilojoule_per_kilogram>()
    );
    // ====================================================================
    // Test (4): Pressure Match
    // ====================================================================
    let p2_result = outlet_state.get_pressure();

    approx::assert_relative_eq!(
        p2.get::<pascal>(),
        p2_result.get::<pascal>(),
        epsilon = 1e-3  // Very tight tolerance for pressure
    );

    println!("✓ Pressure match: p₂ = {:.2} bar", 
        p2_result.get::<bar>()
    );

}

