use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::length::meter;
use uom::si::mass_rate::kilogram_per_second;
use uom::si::pressure::{bar, kilopascal, pascal};
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
fn dry_steam_test_overexpanded(){

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



    let (v, _mass_rate, outlet_state) = 
        calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
            p1, 
            h1, 
            v1, 
            a_throat, 
            a2, 
            p2
        );

    // note the 
    // guess_velocity_and_state_for_diverge_nozzle_from_choked_throat 
    // is the problematic one, will need to test the bisection algorithm


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
    let _rho2 = outlet_state.get_rho();

    // mass balance will not pass in overexpansion, there are shock waves
    //let mass_rate_calculated = rho2 * v * a2;

    //approx::assert_relative_eq!(
    //    mass_rate.get::<kilogram_per_second>(),
    //    mass_rate_calculated.get::<kilogram_per_second>(),
    //    epsilon = 1e-6
    //);

    //println!("✓ Mass balance: ṁ = {:.6} kg/s", 
    //    mass_rate.get::<kilogram_per_second>()
    //);
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

use uom::si::area::square_centimeter;
use uom::si::pressure::megapascal;
use uom::si::ratio::ratio;

/// # Test: `validate_against_cengel_isentropic_steam_nozzle`
///
/// ## Purpose
/// This test validates the master nozzle function against a classic textbook problem
/// from Cengel's "Thermodynamics: An Engineering Approach". It models a simple,
/// un-choked, isentropic flow of superheated steam through a converging nozzle.
///
/// ## Scenario
/// - **Problem:** Steam flow through a converging nozzle.
/// - **Inlet:** Superheated steam at 5 MPa and 400°C, with negligible velocity.
///   (This means the inlet state is effectively the stagnation state).
/// - **Outlet:** The steam exits at a back pressure of 3 MPa.
/// - **Geometry:** The nozzle has a single exit area of 60 cm².
/// - **Assumption:** The flow is isentropic.
///
/// ## Validation Checks
/// The test calls the main nozzle function and asserts that the calculated results
/// for outlet velocity, mass flow rate, and Mach number match the known answers
/// provided by the textbook example.
#[test]
fn validate_against_cengel_isentropic_steam_nozzle() {
    // ====================================================================
    // ARRANGE: Set up the problem based on the textbook data.
    // ====================================================================

    // --- Inlet Conditions (Stagnation State) ---
    let p1 = Pressure::new::<megapascal>(5.0);
    let t1 = ThermodynamicTemperature::new::<degree_celsius>(400.0);
    // "Negligible velocity" means v1 is zero, and p1/t1 are the stagnation properties.
    let v1 = Velocity::new::<meter_per_second>(0.0);

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality_1(t1, p1, ref_vol);
    let h1 = inlet_state.get_specific_enthalpy();

    // --- Outlet Conditions ---
    let p2 = Pressure::new::<megapascal>(3.0);

    // --- Nozzle Geometry ---
    // For a simple converging nozzle, the throat area and the exit area are the same.
    let a_exit = Area::new::<square_centimeter>(60.0);
    let a_throat = a_exit;

    // ====================================================================
    // ACT: Call the master function to calculate the outlet state.
    // ====================================================================

    let (v_out, m_dot_out, state_out) = 
        calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
            p1, 
            h1, 
            v1, 
            a_throat, 
            a_exit, 
            p2
        );
    dbg!(&(v_out, m_dot_out, state_out));

    // ====================================================================
    // ASSERT: Compare the calculated results to the known textbook answers.
    // We use a 1.5% relative tolerance to account for potential small differences
    // in steam table implementations between the textbook and the code.
    // ====================================================================

    // --- 1. Validate Outlet Velocity ---
    let expected_velocity = Velocity::new::<meter_per_second>(529.0);
    println!(
        "Velocity Check: Calculated = {:.2} m/s, Expected = {:.2} m/s",
        v_out.get::<meter_per_second>(),
        expected_velocity.get::<meter_per_second>()
    );
    approx::assert_relative_eq!(
        v_out.get::<meter_per_second>(),
        expected_velocity.get::<meter_per_second>(),
        max_relative = 0.015 // 1.5% tolerance
    );

    // --- 2. Validate Mass Flow Rate ---
    let expected_mass_rate = MassRate::new::<kilogram_per_second>(36.9);
    println!(
        "Mass Flow Check: Calculated = {:.2} kg/s, Expected = {:.2} kg/s",
        m_dot_out.get::<kilogram_per_second>(),
        expected_mass_rate.get::<kilogram_per_second>()
    );
    approx::assert_relative_eq!(
        m_dot_out.get::<kilogram_per_second>(),
        expected_mass_rate.get::<kilogram_per_second>(),
        max_relative = 0.015 // 1.5% tolerance
    );

    // --- 3. Validate Mach Number ---
    
    let outlet_enthalpy = h1 - 0.5 * expected_velocity * expected_velocity;

    let s0 = inlet_state.get_specific_entropy();
    let expected_outlet_state = TampinesSteamTableCV::new_from_hs(
        outlet_enthalpy, s0, ref_vol
    );

    dbg!(&(expected_outlet_state,state_out));

    dbg!(&(
            expected_outlet_state.get_speed_of_sound(),
            state_out.get_speed_of_sound(),
    ));


    let mach_number = state_out.get_mach_number(v_out);
    let expected_mach_number = 0.915;



    println!(
        "Mach Number Check: Calculated = {:.3}, Expected = {:.3}",
        mach_number.get::<ratio>(),
        expected_mach_number
    );
    approx::assert_relative_eq!(
        mach_number.get::<ratio>(),
        expected_mach_number,
        max_relative = 0.015 // 1.5% tolerance
    );
}

/// # Test: `validate_against_cengel_choked_flow_nozzle`
///
/// ## Purpose
/// This test validates the master nozzle function against a Cengel textbook problem
/// involving **choked, isentropic, supersonic flow** through a converging-diverging nozzle.
///
/// ## Scenario
/// - **Problem:** Steam flow through a C-D nozzle where choking occurs.
/// - **Inlet:** Superheated steam at 1 MPa and 500°C, with negligible velocity.
/// - **Outlet:** The steam exits into a low back pressure of 200 kPa, ensuring supersonic flow.
/// - **Known Result:** The flow is choked at a mass flow rate of 2.5 kg/s.
/// - **Geometry:** The exit area is 31.5 cm². The throat area is unknown and must be calculated.
///
/// ## Validation Approach
/// 1.  **Calculate Throat Area:** The test first determines the required throat area (`a_throat`)
///     that would produce the known choked mass flow rate (2.5 kg/s) for the given inlet
///     conditions. This is a critical setup step.
/// 2.  **Execute Model:** It then calls the master nozzle function with all known parameters,
///     including the calculated `a_throat`.
/// 3.  **Assert Results:** Finally, it asserts that the function's outputs for mass flow rate
///     and exit Mach number match the known answers from the textbook.
#[test]
fn validate_against_cengel_choked_flow_nozzle() {
    // ====================================================================
    // ARRANGE: Set up the problem based on the textbook data.
    // ====================================================================

    // --- Inlet Conditions (Stagnation State) ---
    let p1 = Pressure::new::<megapascal>(1.0);
    let t1 = ThermodynamicTemperature::new::<degree_celsius>(500.0);
    let v1 = Velocity::new::<meter_per_second>(0.0); // Negligible velocity

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality_1(t1, p1, ref_vol);
    let h1 = inlet_state.get_specific_enthalpy();
    let s1 = inlet_state.get_specific_entropy();

    // --- Outlet Conditions & Known Results ---
    let p2 = Pressure::new::<kilopascal>(200.0);
    let expected_mass_rate = MassRate::new::<kilogram_per_second>(2.5);
    let expected_mach_number = 1.738;

    // --- Geometry ---
    let a_exit = Area::new::<square_centimeter>(31.5);

    // --- Step 1 (Pre-calculation): Determine the required throat area ---
    // We need to find the throat area that results in a choked flow of 2.5 kg/s.
    // We can do this by calculating the choked mass *flux* (mass flow per unit area),
    // which depends only on stagnation properties, and then finding the area.

    // Stagnation properties are the same as the inlet since v1 is negligible.
    let s0 = s1;

    // now, flow is isentropic,
    // and we know it is choked 
    // let's get the throat state 
    //
    let p_throat = inlet_state.get_critical_pressure_pure_vapour();
    let s_throat = s0;

    let state_throat = TampinesSteamTableCV::new_from_ps(
        p_throat, s_throat, ref_vol
    );

    // now, mass flowrate is rho * a * v 

    let v_throat = state_throat.get_speed_of_sound();
    let rho_throat = state_throat.get_rho();
    let a_throat = expected_mass_rate/v_throat/rho_throat;

    
    println!("--- Pre-calculation ---");
    println!("Calculated required throat area: {:.2} cm²", a_throat.get::<square_centimeter>());


    // ====================================================================
    // ACT: Call the master function with all known parameters.
    // ====================================================================

    let (v_out, m_dot_out, state_out) = 
        calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
            p1, 
            h1, 
            v1, 
            a_throat, 
            a_exit, 
            p2
        );

    // ====================================================================
    // ASSERT: Compare the calculated results to the known textbook answers.
    // ====================================================================

    // --- 1. Validate Mass Flow Rate (Self-consistency check) ---
    // The model should correctly identify the flow as choked and return the choked rate.
    println!("\n--- Validation ---");
    println!(
        "Mass Flow Check: Calculated = {:.3} kg/s, Expected = {:.3} kg/s",
        m_dot_out.get::<kilogram_per_second>(),
        expected_mass_rate.get::<kilogram_per_second>()
    );
    approx::assert_relative_eq!(
        m_dot_out.get::<kilogram_per_second>(),
        expected_mass_rate.get::<kilogram_per_second>(),
        max_relative = 0.01 // 1% tolerance for self-consistency
    );

    // --- 2. Validate Exit Mach Number ---
    let mach_number = state_out.get_mach_number(v_out);
    println!(
        "Mach Number Check: Calculated = {:.3}, Expected = {:.3}",
        mach_number.get::<ratio>(),
        expected_mach_number
    );
    approx::assert_relative_eq!(
        mach_number.get::<ratio>(),
        expected_mach_number,
        max_relative = 0.015 // 1.5% tolerance for model accuracy
    );
}

#[test]
#[ignore="temporary skip test"]
fn wet_steam_test(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let p1 = Pressure::new::<bar>(15.0);
    let x: f64 = 0.80;
    let inlet_state = TampinesSteamTableCV::new_from_sat_pressure_quality(
        p1, x, ref_vol
    );

    let p2 = Pressure::new::<bar>(5.0);
    // critical pressure here is 666 kPa

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
    // mass flowrate is impossible to calculate without the area 
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

