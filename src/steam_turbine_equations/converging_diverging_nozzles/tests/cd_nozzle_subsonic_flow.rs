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

/// # Test: `dry_steam_test`
///
/// ## Purpose
/// This test validates the behavior of the master nozzle function for a simple,
/// **un-choked, subsonic flow** scenario using **superheated (dry) steam**.
///
/// ## Scenario
/// - **Inlet:** Superheated steam at 300°C and 10 bar.
/// - **Outlet:** A high back pressure of 9 bar is set.
/// - **Expected Outcome:** The pressure drop is insufficient to cause choked flow at the
///   throat. The model should correctly identify this as a simple isentropic expansion
///   and calculate the corresponding subsonic outlet state.
///
/// ## Validation Checks
/// The test confirms the physical correctness of the result by verifying four
/// fundamental conservation laws and model assumptions:
///
/// 1.  **Isentropic Flow (`s₁ ≈ s₂`):** Asserts that the specific entropy at the outlet
///     is the same as the inlet, which is the core assumption for this flow regime.
/// 2.  **Mass Balance (`ṁ = ρ₂v₂A₂`):** A self-consistency check confirming that the
///     returned mass flow rate matches the value calculated from the outlet state.
/// 3.  **Energy Balance (`h₀ ≈ h₂ + v₂²/2`):** Asserts that the stagnation enthalpy
///     is conserved throughout the process, confirming the First Law of Thermodynamics.
/// 4.  **Pressure Match (`p₂_actual ≈ p₂_input`):** Verifies that the thermodynamic state
///     returned by the function corresponds to the specified back pressure.
#[test]
fn dry_steam_test(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let temperature = ThermodynamicTemperature::new::<degree_celsius>(300.0);
    let p1 = Pressure::new::<bar>(10.0);
    let inlet_state = TampinesSteamTableCV::new_from_tp_quality_1(
        temperature, p1, ref_vol
    );
    let p2 = Pressure::new::<bar>(9.0);

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
    approx::assert_relative_eq!(
        s2.get::<kilojoule_per_kilogram_kelvin>(),
        s0.get::<kilojoule_per_kilogram_kelvin>(),
        epsilon = 0.001  // 0.1% tolerance for numerical errors
    );
    // For subsonic isentropic flow: s2 = s0 = s1
    approx::assert_relative_eq!(
        s2.get::<kilojoule_per_kilogram_kelvin>(),
        s1.get::<kilojoule_per_kilogram_kelvin>(),
        epsilon = 0.001  // 0.1% tolerance for numerical errors
    );

    println!("✓ Entropy conserved: s1 = {:.4} kJ/kg·K, s2 = {:.4} kJ/kg·K", 
        s1.get::<kilojoule_per_kilogram_kelvin>(),
        s2.get::<kilojoule_per_kilogram_kelvin>()
    );
    
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


/// # Test: `wet_steam_test`
///
/// ## Purpose
/// This test validates the behavior of the master nozzle function for an
/// **un-choked, subsonic flow** scenario, but this time using a **two-phase (wet) steam**
/// mixture.
///
/// ## Scenario
/// - **Inlet:** A saturated steam mixture at 10 bar with a quality of 80%.
/// - **Outlet:** A high back pressure of 9 bar is set.
/// - **Expected Outcome:** Similar to the `dry_steam_test`, the flow should remain
///   subsonic and un-choked. The test verifies that the model handles the isentropic
///   expansion correctly even when starting from a two-phase state.
///
/// ## Validation Checks
/// The test performs the same four fundamental checks as `dry_steam_test` to ensure
/// the model's physical consistency across different fluid phases:
///
/// 1.  **Isentropic Flow (`s₁ ≈ s₂`):** Verifies that entropy is conserved.
/// 2.  **Mass Balance (`ṁ = ρ₂v₂A₂`):** Ensures self-consistency of the returned state.
/// 3.  **Energy Balance (`h₀ ≈ h₂ + v₂²/2`):** Confirms conservation of stagnation enthalpy.
/// 4.  **Pressure Match (`p₂_actual ≈ p₂_input`):** Ensures the outlet state matches the
///     specified back pressure.
#[test]
fn wet_steam_test(){

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let p1 = Pressure::new::<bar>(10.0);
    let x: f64 = 0.80;
    let inlet_state = TampinesSteamTableCV::new_from_sat_pressure_quality(
        p1, x, ref_vol
    );

    let p2 = Pressure::new::<bar>(9.0);

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
    approx::assert_relative_eq!(
        s2.get::<kilojoule_per_kilogram_kelvin>(),
        s0.get::<kilojoule_per_kilogram_kelvin>(),
        epsilon = 0.001  // 0.1% tolerance for numerical errors
    );
    // For subsonic isentropic flow: s2 = s0 = s1
    approx::assert_relative_eq!(
        s2.get::<kilojoule_per_kilogram_kelvin>(),
        s1.get::<kilojoule_per_kilogram_kelvin>(),
        epsilon = 0.001  // 0.1% tolerance for numerical errors
    );

    println!("✓ Entropy conserved: s1 = {:.4} kJ/kg·K, s2 = {:.4} kJ/kg·K", 
        s1.get::<kilojoule_per_kilogram_kelvin>(),
        s2.get::<kilojoule_per_kilogram_kelvin>()
    );
    
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
