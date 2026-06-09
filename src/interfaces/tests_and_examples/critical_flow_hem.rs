use uom::si::area::square_foot;
use uom::si::available_energy::btu_it_per_pound;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::mass_rate::pound_per_second;
use uom::si::pressure::pound_force_per_square_inch;

use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
///// @misc{claude2026tampines,
//  author       = {{Anthropic}},
//  title        = {Conversation on {TAMPINES} {HEM} Critical Flow Implementation},
//  year         = {2026},
//  month        = {June},
//  note         = {AI assistant conversation with Claude Sonnet 4.6 via claude.ai. 
//                  Topics covered: Moody HEM critical mass flux, two-phase 
//                  equilibrium speed of sound, IAPWS-IF97 steam tables in Rust},
//  url          = {https://claude.ai}
//}
/// This is AI generated helper function
/// A reusable helper function to validate any isobar from the Moody chart

/// # Helper Function: `validate_moody_isobar`
///
/// ## Purpose
/// A reusable test function designed to validate the `get_stagnation_critical_mass_flux`
/// method against any given isobar curve from F.J. Moody's 1975 paper.
///
/// ## How it Works
/// 1.  Defines the reference values (`p_ref`, `h_ref`, `g_ref`) used by Moody.
/// 2.  Loops through a provided vector of dimensionless data points `(h/h_ref, G/G_ref)`.
/// 3.  For each point, it reconstructs the physical stagnation state `(p₀, h₀)`.
/// 4.  It calls the model to get the calculated critical mass flux (`g_test`).
/// 5.  It asserts that the calculated flux is within a given tolerance of the
///     theoretical flux from the Moody chart.
fn test_moody_isobar(
    dimensionless_stagnation_pressure: f64,
    data_points: &[(f64, f64)],
    tolerance: f64,
) {
    // --- Define the Reference Values from the Moody Paper ---
    let p_ref = Pressure::new::<pound_force_per_square_inch>(100.0);
    // Note: Moody's paper uses BTU(IT)/lbm, which is what btu_it_per_pound represents.
    let h_ref = AvailableEnergy::new::<btu_it_per_pound>(100.0);
    let g_ref: MassFlux = MassRate::new::<pound_per_second>(1000.0) / Area::new::<square_foot>(1.0);
    let ref_vol = TampinesSteamTableCV::get_ref_vol();

    // --- Loop Through Each Data Point for the Given Isobar ---
    for (h_dimensionless_ptr, g_dimensionless_ptr) in data_points.iter() {
        let h0 = h_ref * (*h_dimensionless_ptr);
        let p0 = dimensionless_stagnation_pressure * p_ref;
        let g_ref_expected = g_ref * (*g_dimensionless_ptr);

        let state_0 = TampinesSteamTableCV::new_from_ph(p0, h0, ref_vol);
        let g_test = state_0.get_stagnation_critical_mass_flux();
        // this helps see which point we are at on the graph
        dbg!(&(*h_dimensionless_ptr,*g_dimensionless_ptr,g_test/g_ref));
        dbg!(&(g_ref_expected,g_test));

        // The assertion uses the provided tolerance to compare the model's result
        // against the theoretical value from the Moody chart.
        approx::assert_relative_eq!(
            g_ref_expected.get::<kilogram_per_square_meter_second>().log10(),
            g_test.get::<kilogram_per_square_meter_second>().log10(),
            max_relative = tolerance
        );
    }
}
/// # Test: `isobar_pref_12_00`
/// Validates the critical mass flux model against the `p/p_ref = 12.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_12_00() {
    let data = vec![
        //(0.6471, 26.6543), 
        //(1.0588, 26.0586), 
        //(1.5098, 26.0586), 
        //(2.0784, 25.4762),
        //(2.5882, 24.627), 
        //(2.9608, 23.806), 
        //(3.5098, 23.0125), 
        //(3.9216, 21.0232),
        //(4.3137, 19.2059), 
        //(4.6667, 17.1536), 
        //(5.0196, 14.3162), (5.1961, 11.8139),
        (5.3529, 9.1098), (5.5294, 6.7905), (5.8039, 5.6036), (6.2353, 4.7836),
        (6.8235, 4.13), (7.4118, 3.6887), (8.0784, 3.3698), (8.6863, 3.1135),
        (9.3333, 2.8444), (9.8431, 2.7187), (10.4314, 2.5985), (11.098, 2.4557),
        (11.7451, 2.3739),
    ];
    test_moody_isobar(12.00, &data, 1e-2);
}

