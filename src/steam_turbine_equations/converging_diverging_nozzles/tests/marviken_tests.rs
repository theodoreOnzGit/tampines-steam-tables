// for the Marviken tests 
//
// there was a vessel at 425 m3  (fig 3.1 on page 9)
// https://www.nrc.gov/docs/ML2005/ML20052H367.pdf
// This vessel was 

// In your test module

use uom::si::f64::*;
use uom::si::mass_rate::kilogram_per_second;
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::pressure::{atmosphere, megapascal};
use uom::si::length::{meter, millimeter};
use uom::si::velocity::meter_per_second;

use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
use crate::steam_turbine_equations::calculate_velocity_mass_flowrate_and_state_in_cd_nozzle;

/// looks like fig 8:24 on page 113 seems to be the best 
/// as most nozzles have pipe which may cause extra pressure loss and 
/// reduce flowrate.
///
/// These correspond to test 23 and 24. 
///
/// But Fig 8:24 has L/D at 0.3 
/// which makes it as pure a nozzle as it gets
///
/// In the pdf, this is on page 100 (113 of the pdf)
/// The conditions of test 23 and 24 are in table 4:2 on 
/// page 21 (page 34 of the pdf)
///
/// https://www.nrc.gov/docs/ML2005/ML20052H367.pdf
/// NUREG/CR-2671
/// MXC-301
///
/// Theri
#[test]
fn validate_against_marviken_test_24() {
    // --- Step 1: Define Initial Conditions from Table 4.1.1 ---
    let p1 = Pressure::new::<megapascal>(4.95);
    let t1 = TampinesSteamTableCV::try_get_tsat(p1).unwrap() - 
        TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(33.0);
    
    // The back pressure is atmospheric, as they are venting to a large containment vessel.
    let p2 = Pressure::new::<atmosphere>(1.0);

    // --- Step 2: Define Geometry from Figure 3.2.2 ---
    let nozzle_diameter = Length::new::<millimeter>(500.0);
    let nozzle_area = std::f64::consts::PI * (nozzle_diameter * nozzle_diameter / 4.0);
    let a_throat = nozzle_area;
    let a_exit = nozzle_area; // It's a converging nozzle, so throat area = exit area

    // --- Step 3: Get Initial State from your Steam Tables ---
    let ref_vol = TampinesSteamTableCV::get_ref_vol();
    let state_1 = TampinesSteamTableCV::new_from_tp_quality_1(t1, p1, ref_vol);
    let h1 = state_1.get_specific_enthalpy();

    // The velocity inside the huge Marviken pressure vessel is effectively zero.
    let v1 = Velocity::new::<meter_per_second>(0.0);

    // --- Step 4: Call Your Master Function ---
    let (v_out, m_dot_out, state_out) = 
        calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
            p1, h1, v1, a_throat, a_exit, p2
        );

    // --- Step 5: Compare to the Experimental Result from Table 5.1.1 ---
    let experimental_mass_flowrate = MassRate::new::<kilogram_per_second>(895.0);

    println!("Calculated Mass Flow Rate: {:.2} kg/s", m_dot_out.get::<kilogram_per_second>());
    println!("Experimental Mass Flow Rate: {:.2} kg/s", experimental_mass_flowrate.get::<kilogram_per_second>());

    // Use a relative difference assertion. Don't expect a perfect match!
    // Getting within 10-15% would be a fantastic result for a 1D model.
    let relative_difference = ((m_dot_out - experimental_mass_flowrate) / experimental_mass_flowrate).abs();
    
    println!("Relative Difference: {:.2}%", relative_difference.get::<ratio>() * 100.0);

    // Assert that your result is within a reasonable tolerance, e.g., 20%
    assert!(relative_difference.get::<ratio>() < 0.20);
}
