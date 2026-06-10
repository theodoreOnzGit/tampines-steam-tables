use uom::si::area::square_foot;
use uom::si::available_energy::btu_it_per_pound;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::mass_rate::pound_per_second;
use uom::si::pressure::pound_force_per_square_inch;

use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
use crate::steam_turbine_equations::choked_flow::{g_max_hem_analytical_ph, isentropic_pressure_scan_of_mass_flux};

// please note for the test:
// For p0/p_ref = 0.25
// p0 = 1.724 bar (0.78% of p_crit = 220.64 bar)
// T_sat ≈ 115.2°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 0.50
// p0 = 3.447 bar (1.56% of p_crit)
// T_sat ≈ 138.9°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 1.00
// p0 = 6.895 bar (3.13% of p_crit)
// T_sat ≈ 165.0°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 2.00
// p0 = 13.790 bar (6.25% of p_crit)
// T_sat ≈ 191.6°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 4.00
// p0 = 27.579 bar (12.50% of p_crit)
// T_sat ≈ 230.1°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 6.00
// p0 = 41.369 bar (18.75% of p_crit)
// T_sat ≈ 253.3°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 8.00
// p0 = 55.158 bar (25.00% of p_crit)
// T_sat ≈ 270.3°C
// Inlet region: Region 4 (wet steam) or Region 2 (superheated steam)
// #[ignore] - failing

// For p0/p_ref = 10.00
// p0 = 68.948 bar (31.25% of p_crit)
// T_sat ≈ 284.6°C
// Inlet region: Region 4 (wet steam) or Region 1 (subcooled liquid)
// #[ignore] - failing

// For p0/p_ref = 12.00
// p0 = 82.737 bar (37.50% of p_crit)
// T_sat ≈ 296.9°C
// Inlet region: Region 4 (wet steam) or Region 1 (subcooled liquid)
// #[ignore] - failing

// For p0/p_ref = 14.00
// p0 = 96.527 bar (43.75% of p_crit)
// T_sat ≈ 307.8°C
// Inlet region: Region 1 (subcooled liquid) for low h, Region 4 for higher h
// passing

// For p0/p_ref = 16.00
// p0 = 110.316 bar (50.00% of p_crit)
// T_sat ≈ 317.6°C
// Inlet region: Region 1 (subcooled liquid) for low h, Region 4 for higher h
// passing

// For p0/p_ref = 20.00
// p0 = 137.895 bar (62.50% of p_crit)
// T_sat ≈ 335.2°C
// Inlet region: Region 1 (subcooled liquid) for low h, Region 4 for higher h
// passing

// For p0/p_ref = 30.00
// p0 = 206.843 bar (93.74% of p_crit)
// T_sat ≈ 365.8°C — near critical point, Region 3 behaviour
// Inlet region: Region 1 (subcooled liquid) for low h, Region 3 near critical
// passing

/// From Figure 1 of:
///
/// Moody, F. J. (1975). Maximum discharge rate of liquid-vapor mixtures 
/// from vessels (No.
/// NEDO--21052). General Electric Co., San Jose, CA (United States). 
/// BWR Projects Dept..0 
///
/// Downloaded at:
/// https://www.osti.gov/servlets/purl/7309475
///
/// Based on Kretzchmar,
/// we can use eq 2.80 and table 2.140 to calculate any property 
/// desired directly rather than through iteration 
/// Unsure whether this is helpful, but it is a clue...
///
///
/// p0/p_ref = 0.25
///"dimensionless stagnation enthalpy","dimensionless critical mass flux"
/// 0.4902,3.8593
/// 0.8039,3.7306
/// 1.1961,3.4469
/// 1.4314,3.1135
/// 1.6275,2.7187
/// 1.7647,2.2948
/// 1.8627,1.5106
/// 1.902,1.1133
/// 1.9412,0.4991
/// 1.9804,0.3678
/// 2.0588,0.2901
/// 2.2157,0.2212
/// 2.4118,0.1867
/// 2.8431,0.1541
/// 3.1176,0.1392
/// 3.5098,0.1243
/// 3.8431,0.111
/// 4.4902,0.1014
/// 5.0784,0.0916
/// 5.5098,0.0866
/// 6.1765,0.0791
/// 6.6863,0.0756
/// 7.2549,0.0706
/// 7.8039,0.0691
/// 8.3725,0.0653
/// 8.7255,0.0631
/// 9.3333,0.061
/// 9.902,0.057
/// 10.1765,0.057
/// 10.549,0.0564
/// 11.0784,0.0551
/// 11.5098,0.0533
/// 
///
#[test]
#[ignore]
fn isobar_pref_0_25() {


    // let's first have our datapoints 

    let isobar_pref_0_25_critical_dimensionless_critical_mass_flux = vec![
        (0.4902,3.8593),
        (0.8039,3.7306),
        (1.1961,3.4469),
        (1.4314,3.1135),
        (1.6275,2.7187),
        (1.7647,2.2948),
        (1.8627,1.5106),
        (1.902,1.1133),
        (1.9412,0.4991),
        (1.9804,0.3678),
        (2.0588,0.2901),
        (2.2157,0.2212),
        (2.4118,0.1867),
        (2.8431,0.1541),
        (3.1176,0.1392),
        (3.5098,0.1243),
        (3.8431,0.111),
        (4.4902,0.1014),
        (5.0784,0.0916),
        (5.5098,0.0866),
        (6.1765,0.0791),
        (6.6863,0.0756),
        (7.2549,0.0706),
        (7.8039,0.0691),
        (8.3725,0.0653),
        (8.7255,0.0631),
        (9.3333,0.061),
        (9.902,0.057),
        (10.1765,0.057),
        (10.549,0.0564),
        (11.0784,0.0551),
        (11.5098,0.0533),
        ];


    let p_ref = Pressure::new::<pound_force_per_square_inch>(100.0);
    let dimensionless_stagnation_pressure = 0.25;
    // for clarity: it is 2.326e5 joule/kg
    // it matches closer to btu_it as opposed to 
    // btu, which is:
    // 2.324_443_707_610_621_E3 joule/kg
    let h_ref = AvailableEnergy::new::<btu_it_per_pound>(100.0);
    let g_ref: MassFlux = 
        MassRate::new::<pound_per_second>(1000.0)/
        Area::new::<square_foot>(1.0);

    // ref vol 
    let ref_vol = TampinesSteamTableCV::get_ref_vol();

    for (h_dimensionless_ptr, g_dimensionless_ptr) 
        in isobar_pref_0_25_critical_dimensionless_critical_mass_flux.iter() {

            let h0: AvailableEnergy = h_ref * (*h_dimensionless_ptr);
            let p0: Pressure = dimensionless_stagnation_pressure * p_ref;

            let g_ref_expected: MassFlux = g_ref * (*g_dimensionless_ptr);

            let state_0 = TampinesSteamTableCV::new_from_ph(
                p0, h0, ref_vol
            );

            let g_test = state_0.get_stagnation_critical_mass_flux();
            // this helps see which point we are at on the graph
            dbg!(&(*h_dimensionless_ptr,*g_dimensionless_ptr,g_test/g_ref));

            // note: I took these values from a log (y) vs x graph 
            // as in log (g_dimensionless) vs h_dimensionless graph 
            // hence, errors will be big on for the larger values, for 
            // graphreader
            // it is better to assert errors on the log scale rather than 
            // the linear scale, until such time I get data from linear 
            // scale graph

            approx::assert_relative_eq!(
                g_ref_expected.get::<kilogram_per_square_meter_second>().log10(),
                g_test.get::<kilogram_per_square_meter_second>().log10(),
                max_relative=1e-2
            );


    }



}

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
fn validate_moody_isobar(
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

        // The assertion uses the provided tolerance to compare the model's result
        // against the theoretical value from the Moody chart.
        approx::assert_relative_eq!(
            g_ref_expected.get::<kilogram_per_square_meter_second>().log10(),
            g_test.get::<kilogram_per_square_meter_second>().log10(),
            max_relative = tolerance
        );
    }
}

///
fn validate_moody_isobar_hem(
    dimensionless_stagnation_pressure: f64,
    data_points: &[(f64, f64)],
    tolerance: f64,
) {
    // --- Define the Reference Values from the Moody Paper ---
    let p_ref = Pressure::new::<pound_force_per_square_inch>(100.0);
    // Note: Moody's paper uses BTU(IT)/lbm, which is what btu_it_per_pound represents.
    let h_ref = AvailableEnergy::new::<btu_it_per_pound>(100.0);
    let g_ref: MassFlux = MassRate::new::<pound_per_second>(1000.0) / Area::new::<square_foot>(1.0);

    // --- Loop Through Each Data Point for the Given Isobar ---
    for (h_dimensionless_ptr, g_dimensionless_ptr) in data_points.iter() {
        let h0 = h_ref * (*h_dimensionless_ptr);
        let p0 = dimensionless_stagnation_pressure * p_ref;
        let g_ref_expected = g_ref * (*g_dimensionless_ptr);

        let g_test = g_max_hem_analytical_ph(p0, h0);
        // this helps see which point we are at on the graph
        dbg!(&(*h_dimensionless_ptr,*g_dimensionless_ptr,g_test/g_ref));

        // The assertion uses the provided tolerance to compare the model's result
        // against the theoretical value from the Moody chart.
        approx::assert_relative_eq!(
            g_ref_expected.get::<kilogram_per_square_meter_second>().log10(),
            g_test.get::<kilogram_per_square_meter_second>().log10(),
            max_relative = tolerance
        );
    }
}

// For p0/p_ref = 0.50 
// "x","y"
// 0.4902,5.4168
// 0.7647,5.2362
// 1.2353,5.0617
// 1.6471,4.6241
// 1.9412,3.9031
// 2.1765,3.1135
// 2.2549,2.269
// 2.3137,1.275
// 2.4118,0.7005
// 2.6078,0.4508
// 3,0.336
// 3.4314,0.2773
// 3.9804,0.2314
// 4.5686,0.2021
// 5.1373,0.1867
// 5.8235,0.1668
// 6.2157,0.1612
// 7.1569,0.144
// 8.1373,0.133
// 8.8627,0.1271
// 9.8431,0.1148
// 10.5686,0.111
// 11.2549,0.1073
// 11.7255,0.1037


/// # Test: `isobar_pref_0_50`
/// Validates the critical mass flux model against the `p/p_ref = 0.50` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_0_50() {
    let data = vec![
        (0.4902, 5.4168), (0.7647, 5.2362), (1.2353, 5.0617), (1.6471, 4.6241),
        (1.9412, 3.9031), (2.1765, 3.1135), (2.2549, 2.269), (2.3137, 1.275),
        (2.4118, 0.7005), (2.6078, 0.4508), (3.0, 0.336), (3.4314, 0.2773),
        (3.9804, 0.2314), (4.5686, 0.2021), (5.1373, 0.1867), (5.8235, 0.1668),
        (6.2157, 0.1612), (7.1569, 0.144), (8.1373, 0.133), (8.8627, 0.1271),
        (9.8431, 0.1148), (10.5686, 0.111), (11.2549, 0.1073), (11.7255, 0.1037),
    ];
    validate_moody_isobar(0.50, &data, 1e-2);
}
#[test]
#[ignore]
fn isobar_pref_0_50_hem() {
    let data = vec![
        (0.4902, 5.4168), (0.7647, 5.2362), (1.2353, 5.0617), (1.6471, 4.6241),
        (1.9412, 3.9031), (2.1765, 3.1135), (2.2549, 2.269), (2.3137, 1.275),
        (2.4118, 0.7005), (2.6078, 0.4508), (3.0, 0.336), (3.4314, 0.2773),
        (3.9804, 0.2314), (4.5686, 0.2021), (5.1373, 0.1867), (5.8235, 0.1668),
        (6.2157, 0.1612), (7.1569, 0.144), (8.1373, 0.133), (8.8627, 0.1271),
        (9.8431, 0.1148), (10.5686, 0.111), (11.2549, 0.1073), (11.7255, 0.1037),
    ];
    validate_moody_isobar_hem(0.50, &data, 1e-2);
}
#[test]
#[ignore]
fn isobar_pref_0_50_pressure_scan() {
    let data = vec![
        (0.4902, 5.4168), 
        //(0.7647, 5.2362), (1.2353, 5.0617), (1.6471, 4.6241),
        //(1.9412, 3.9031), (2.1765, 3.1135), (2.2549, 2.269), (2.3137, 1.275),
        //(2.4118, 0.7005), (2.6078, 0.4508), (3.0, 0.336), (3.4314, 0.2773),
        //(3.9804, 0.2314), (4.5686, 0.2021), (5.1373, 0.1867), (5.8235, 0.1668),
        //(6.2157, 0.1612), (7.1569, 0.144), (8.1373, 0.133), (8.8627, 0.1271),
        //(9.8431, 0.1148), (10.5686, 0.111), (11.2549, 0.1073), (11.7255, 0.1037),
    ];
    // --- Define the Reference Values from the Moody Paper ---
    let p_ref = Pressure::new::<pound_force_per_square_inch>(100.0);
    let dimensionless_stagnation_pressure = 0.50;
    // Note: Moody's paper uses BTU(IT)/lbm, which is what btu_it_per_pound represents.
    let h_ref = AvailableEnergy::new::<btu_it_per_pound>(100.0);
    let g_ref: MassFlux = MassRate::new::<pound_per_second>(1000.0) / Area::new::<square_foot>(1.0);

    // --- Loop Through Each Data Point for the Given Isobar ---
    for (h_dimensionless_ptr, g_dimensionless_ptr) in data.iter() {
        let h0 = h_ref * (*h_dimensionless_ptr);
        let p0 = dimensionless_stagnation_pressure * p_ref;
        let g_ref_expected = g_ref * (*g_dimensionless_ptr);

        let s0 = s_ph_eqm(p0, h0);

        isentropic_pressure_scan_of_mass_flux(s0, p0);

        // this helps see which point we are at on the graph
        dbg!(&(*h_dimensionless_ptr,g_ref_expected));

        // The assertion uses the provided tolerance to compare the model's result
        // against the theoretical value from the Moody chart.

    }
    todo!();
}

// For p0/p_ref = 1.00
// "x","y"
// 0.451,7.6029
// 0.6667,7.6029
// 1.3137,7.1852
// 1.9412,6.5641
// 2.4314,5.0617
// 2.6471,3.9475
// 2.8235,2.7496
// 2.8627,1.9812
// 2.9216,1.205
// 3.0588,0.8393
// 3.5686,0.598
// 4.2941,0.4559
// 5.1373,0.3762
// 6.0392,0.3248
// 6.6275,0.3001
// 7.451,0.2742
// 8.2549,0.2591
// 9.2353,0.2394
// 10.1373,0.2263
// 11.1961,0.2138
// 11.8235,0.2067


/// # Test: `isobar_pref_1_00`
/// Validates the critical mass flux model against the `p/p_ref = 1.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_1_00() {
    let data = vec![
        (0.451, 7.6029), (0.6667, 7.6029), (1.3137, 7.1852), (1.9412, 6.5641),
        (2.4314, 5.0617), (2.6471, 3.9475), (2.8235, 2.7496), (2.8627, 1.9812),
        (2.9216, 1.205), (3.0588, 0.8393), (3.5686, 0.598), (4.2941, 0.4559),
        (5.1373, 0.3762), (6.0392, 0.3248), (6.6275, 0.3001), (7.451, 0.2742),
        (8.2549, 0.2591), (9.2353, 0.2394), (10.1373, 0.2263), (11.1961, 0.2138),
        (11.8235, 0.2067),
    ];
    validate_moody_isobar(1.00, &data, 1e-2);
}

// For p0/p_ref = 2.0
// "x","y"
// 0.4706,11.0394
// 0.9216,10.7927
// 1.3529,10.4329
// 1.9412,10.1997
// 2.5294,9.0074
// 2.8824,7.2669
// 3.1569,5.7968
// 3.3333,3.9923
// 3.4314,2.4009
// 3.4902,1.6723
// 3.9216,1.1781
// 4.5098,0.9397
// 5.2353,0.7843
// 6.0784,0.6695
// 7.3725,0.5716
// 8.9804,0.4879
// 10.4314,0.4358
// 11.902,0.3981


/// # Test: `isobar_pref_2_00`
/// Validates the critical mass flux model against the `p/p_ref = 2.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_2_00() {
    let data = vec![
        (0.4706, 11.0394), (0.9216, 10.7927), (1.3529, 10.4329), (1.9412, 10.1997),
        (2.5294, 9.0074), (2.8824, 7.2669), (3.1569, 5.7968), (3.3333, 3.9923),
        (3.4314, 2.4009), (3.4902, 1.6723), (3.9216, 1.1781), (4.5098, 0.9397),
        (5.2353, 0.7843), (6.0784, 0.6695), (7.3725, 0.5716), (8.9804, 0.4879),
        (10.4314, 0.4358), (11.902, 0.3981),
    ];
    validate_moody_isobar(2.00, &data, 1e-2);
}

// For p0/p_ref = 4.0 
// "x","y"
// 0.7255,13.2273
// 1.1961,12.7864
// 1.6078,12.7864
// 2.0588,12.3602
// 2.5294,11.9481
// 2.8039,11.1648
// 3.1176,9.6394
// 3.3529,8.4169
// 3.5882,6.5641
// 3.7255,4.7298
// 3.7843,3.332
// 3.902,2.269
// 4.2745,1.7105
// 4.6471,1.4602
// 5.1373,1.2325
// 5.7255,1.1008
// 6.451,0.9832
// 7.451,0.8393
// 8.4118,0.7581
// 9.1569,0.7005
// 10.098,0.6771
// 10.9804,0.6327
// 11.8235,0.6256


/// # Test: `isobar_pref_4_00`
/// Validates the critical mass flux model against the `p/p_ref = 4.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_4_00() {
    let data = vec![
        (0.7255, 13.2273), (1.1961, 12.7864), (1.6078, 12.7864), (2.0588, 12.3602),
        (2.5294, 11.9481), (2.8039, 11.1648), (3.1176, 9.6394), (3.3529, 8.4169),
        (3.5882, 6.5641), (3.7255, 4.7298), (3.7843, 3.332), (3.902, 2.269),
        (4.2745, 1.7105), (4.6471, 1.4602), (5.1373, 1.2325), (5.7255, 1.1008),
        (6.451, 0.9832), (7.451, 0.8393), (8.4118, 0.7581), (9.1569, 0.7005),
        (10.098, 0.6771), (10.9804, 0.6327), (11.8235, 0.6256),
    ];
    validate_moody_isobar(4.00, &data, 1e-2);
}

// For p0/p_ref = 6.0 
//
// "x","y"
// 0.5882,18.5657
// 0.8235,18.5657
// 1.1765,18.5657
// 1.5098,18.3571
// 1.9608,17.9468
// 2.2941,17.7452
// 2.6275,16.9609
// 3.0196,16.3955
// 3.5882,14.1553
// 4.1373,10.6714
// 4.3922,8.0449
// 4.5098,5.7317
// 4.5686,4.321
// 4.7647,3.6062
// 5.098,2.9759
// 5.4902,2.5693
// 5.9216,2.269
// 6.6078,2.0037
// 7.2353,1.8305
// 7.6863,1.7105
// 8.4314,1.5451
// 9.0196,1.4768
// 9.7647,1.3645
// 10.549,1.275
// 11.2549,1.2187
// 11.7843,1.1517


/// # Test: `isobar_pref_6_00`
/// Validates the critical mass flux model against the `p/p_ref = 6.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_6_00() {
    let data = vec![
        (0.5882, 18.5657), (0.8235, 18.5657), (1.1765, 18.5657), (1.5098, 18.3571),
        (1.9608, 17.9468), (2.2941, 17.7452), (2.6275, 16.9609), (3.0196, 16.3955),
        (3.5882, 14.1553), (4.1373, 10.6714), (4.3922, 8.0449), (4.5098, 5.7317),
        (4.5686, 4.321), (4.7647, 3.6062), (5.098, 2.9759), (5.4902, 2.5693),
        (5.9216, 2.269), (6.6078, 2.0037), (7.2353, 1.8305), (7.6863, 1.7105),
        (8.4314, 1.5451), (9.0196, 1.4768), (9.7647, 1.3645), (10.549, 1.275),
        (11.2549, 1.2187), (11.7843, 1.1517),
    ];
    validate_moody_isobar(6.00, &data, 1e-2);
}

// For p0/p_ref = 8.0
// "x","y"
// 0.6471,21.9954
// 0.8824,21.7482
// 1.451,21.2622
// 2.0196,20.787
// 2.5098,20.0941
// 3,19.2059
// 3.451,17.9468
// 3.8824,16.0291
// 4.2549,13.9963
// 4.4706,11.8139
// 4.6863,9.424
// 4.8039,7.6029
// 4.902,5.8627
// 5.0784,4.7298
// 5.3725,3.9031
// 5.7451,3.486
// 6.1961,3.0439
// 6.5686,2.7808
// 7.0196,2.5404
// 7.5294,2.3739
// 7.9412,2.2435
// 8.2745,2.1202
// 8.8039,2.0495
// 9.2941,1.9152
// 9.7059,1.8513
// 10.1765,1.7496
// 10.6471,1.73
// 11.0784,1.6165
// 11.6078,1.5804
// 11.9608,1.5804



/// # Test: `isobar_pref_8_00`
/// Validates the critical mass flux model against the `p/p_ref = 8.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_8_00() {
    let data = vec![
        (0.6471, 21.9954), (0.8824, 21.7482), (1.451, 21.2622), (2.0196, 20.787),
        (2.5098, 20.0941), (3.0, 19.2059), (3.451, 17.9468), (3.8824, 16.0291),
        (4.2549, 13.9963), (4.4706, 11.8139), (4.6863, 9.424), (4.8039, 7.6029),
        (4.902, 5.8627), (5.0784, 4.7298), (5.3725, 3.9031), (5.7451, 3.486),
        (6.1961, 3.0439), (6.5686, 2.7808), (7.0196, 2.5404), (7.5294, 2.3739),
        (7.9412, 2.2435), (8.2745, 2.1202), (8.8039, 2.0495), (9.2941, 1.9152),
        (9.7059, 1.8513), (10.1765, 1.7496), (10.6471, 1.73), (11.0784, 1.6165),
        (11.6078, 1.5804), (11.9608, 1.5804),
    ];
    validate_moody_isobar(8.00, &data, 1e-2);
}

// For p0/p_ref = 10.0 
// "x","y"
// 0.6275,24.627
// 0.902,24.0766
// 1.2941,24.3502
// 1.8627,23.5385
// 2.2353,23.0125
// 3.0392,21.7482
// 3.5882,20.0941
// 4.1569,17.9468
// 4.5686,15.3206
// 4.8627,12.6427
// 5.0196,9.9718
// 5.1176,8.1363
// 5.2549,6.2035
// 5.5686,4.9486
// 6.0392,4.1769
// 6.8235,3.4469
// 7.5294,3.0439
// 8.2745,2.6881
// 9.2353,2.4282
// 10,2.2183
// 10.8431,2.0964
// 11.4118,2.0495
// 11.7451,2.0037


/// # Test: `isobar_pref_10_00`
/// Validates the critical mass flux model against the `p/p_ref = 10.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_10_00() {
    let data = vec![
        (0.6275, 24.627), (0.902, 24.0766), (1.2941, 24.3502), (1.8627, 23.5385),
        (2.2353, 23.0125), (3.0392, 21.7482), (3.5882, 20.0941), (4.1569, 17.9468),
        (4.5686, 15.3206), (4.8627, 12.6427), (5.0196, 9.9718), (5.1176, 8.1363),
        (5.2549, 6.2035), (5.5686, 4.9486), (6.0392, 4.1769), (6.8235, 3.4469),
        (7.5294, 3.0439), (8.2745, 2.6881), (9.2353, 2.4282), (10.0, 2.2183),
        (10.8431, 2.0964), (11.4118, 2.0495), (11.7451, 2.0037),
    ];
    validate_moody_isobar(10.00, &data, 1e-2);
}

// For p0/p_ref = 12.0 
// "x","y"
// 0.6471,26.6543
// 1.0588,26.0586
// 1.5098,26.0586
// 2.0784,25.4762
// 2.5882,24.627
// 2.9608,23.806
// 3.5098,23.0125
// 3.9216,21.0232
// 4.3137,19.2059
// 4.6667,17.1536
// 5.0196,14.3162
// 5.1961,11.8139
// 5.3529,9.1098
// 5.5294,6.7905
// 5.8039,5.6036
// 6.2353,4.7836
// 6.8235,4.13
// 7.4118,3.6887
// 8.0784,3.3698
// 8.6863,3.1135
// 9.3333,2.8444
// 9.8431,2.7187
// 10.4314,2.5985
// 11.098,2.4557
// 11.7451,2.3739


/// # Test: `isobar_pref_12_00`
/// Validates the critical mass flux model against the `p/p_ref = 12.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_12_00() {
    let data = vec![
        (0.6471, 26.6543), (1.0588, 26.0586), (1.5098, 26.0586), (2.0784, 25.4762),
        (2.5882, 24.627), (2.9608, 23.806), (3.5098, 23.0125), (3.9216, 21.0232),
        (4.3137, 19.2059), (4.6667, 17.1536), (5.0196, 14.3162), (5.1961, 11.8139),
        (5.3529, 9.1098), (5.5294, 6.7905), (5.8039, 5.6036), (6.2353, 4.7836),
        (6.8235, 4.13), (7.4118, 3.6887), (8.0784, 3.3698), (8.6863, 3.1135),
        (9.3333, 2.8444), (9.8431, 2.7187), (10.4314, 2.5985), (11.098, 2.4557),
        (11.7451, 2.3739),
    ];
    validate_moody_isobar(12.00, &data, 1e-2);
}

// For p0/p_ref = 14.0 
// "x","y"
// 0.6667,28.8485
// 1.0588,28.2037
// 1.5294,28.2037
// 1.9412,28.2037
// 2.4706,26.9572
// 3.098,25.7658
// 3.6471,24.0766
// 4.2353,22.2454
// 4.7059,19.645
// 5.0588,16.9609
// 5.4118,13.0787
// 5.6078,10.0851
// 5.7255,8.3223
// 5.9608,6.5641
// 6.2353,5.7968
// 6.5686,5.3559
// 7,4.8929
// 7.4706,4.5208
// 8.0392,4.0377
// 8.549,3.773
// 9.2941,3.3698
// 9.9804,3.1135
// 10.6078,2.9425
// 11.2549,2.8444
// 11.7255,2.7496


/// # Test: `isobar_pref_14_00`
/// Validates the critical mass flux model against the `p/p_ref = 14.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_14_00() {
    let data = vec![
        (0.6667, 28.8485), (1.0588, 28.2037), (1.5294, 28.2037), (1.9412, 28.2037),
        (2.4706, 26.9572), (3.098, 25.7658), (3.6471, 24.0766), (4.2353, 22.2454),
        (4.7059, 19.645), (5.0588, 16.9609), (5.4118, 13.0787), (5.6078, 10.0851),
        (5.7255, 8.3223), (5.9608, 6.5641), (6.2353, 5.7968), (6.5686, 5.3559),
        (7.0, 4.8929), (7.4706, 4.5208), (8.0392, 4.0377), (8.549, 3.773),
        (9.2941, 3.3698), (9.9804, 3.1135), (10.6078, 2.9425), (11.2549, 2.8444),
        (11.7255, 2.7496),
    ];
    validate_moody_isobar(14.00, &data, 1e-2);
}

// For p0/p_ref = 16.0 
// "x","y"
// 0.7059,30.1825
// 1.0196,30.1825
// 1.3333,30.1825
// 1.7255,30.1825
// 2.098,29.5079
// 2.5098,29.1763
// 2.8431,28.5243
// 3.2549,27.5734
// 3.6471,26.3547
// 4.0196,25.1899
// 4.4314,23.0125
// 4.8235,21.2622
// 5.1373,19.2059
// 5.3922,16.3955
// 5.6471,13.6835
// 5.8039,11.0394
// 6,8.7072
// 6.3137,6.7905
// 6.6863,6.0649
// 7.1961,5.5406
// 7.7059,4.9486
// 8.3137,4.5208
// 8.902,4.1769
// 9.5294,3.773
// 10.3137,3.5257
// 10.9804,3.3698
// 11.5882,3.2209


/// # Test: `isobar_pref_16_00`
/// Validates the critical mass flux model against the `p/p_ref = 16.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_16_00() {
    let data = vec![
        (0.7059, 30.1825), (1.0196, 30.1825), (1.3333, 30.1825), (1.7255, 30.1825),
        (2.098, 29.5079), (2.5098, 29.1763), (2.8431, 28.5243), (3.2549, 27.5734),
        (3.6471, 26.3547), (4.0196, 25.1899), (4.4314, 23.0125), (4.8235, 21.2622),
        (5.1373, 19.2059), (5.3922, 16.3955), (5.6471, 13.6835), (5.8039, 11.0394),
        (6.0, 8.7072), (6.3137, 6.7905), (6.6863, 6.0649), (7.1961, 5.5406),
        (7.7059, 4.9486), (8.3137, 4.5208), (8.902, 4.1769), (9.5294, 3.773),
        (10.3137, 3.5257), (10.9804, 3.3698), (11.5882, 3.2209),
    ];
    validate_moody_isobar(16.00, &data, 1e-2);
}

// For p0/p_ref = 20.0 
// "x","y"
// 0.6863,34.1777
// 0.9608,34.1777
// 1.4118,34.1777
// 1.8235,33.7936
// 2.2353,33.4138
// 2.6863,32.6671
// 3.098,31.2233
// 3.6667,29.8433
// 4.0784,28.5243
// 4.4902,26.6543
// 4.8824,24.627
// 5.1569,22.4982
// 5.4902,20.0941
// 5.8039,17.5457
// 6.0392,14.8099
// 6.2745,11.9481
// 6.3922,9.8597
// 6.7059,8.1363
// 6.9608,7.5175
// 7.3333,6.9457
// 7.8431,6.274
// 8.3137,5.8627
// 8.8235,5.4168
// 9.3333,5.0617
// 9.9216,4.6767
// 10.3725,4.47
// 10.9608,4.2244
// 11.2941,4.13
/// # Test: `isobar_pref_20_00`
/// Validates the critical mass flux model against the `p/p_ref = 20.00` isobar
/// from Figure 1 of Moody (1975).
#[test]
#[ignore]
fn isobar_pref_20_00() {
    let data = vec![
        (0.6863, 34.1777), (0.9608, 34.1777), (1.4118, 34.1777), (1.8235, 33.7936),
        (2.2353, 33.4138), (2.6863, 32.6671), (3.098, 31.2233), (3.6667, 29.8433),
        (4.0784, 28.5243), (4.4902, 26.6543), (4.8824, 24.627), (5.1569, 22.4982),
        (5.4902, 20.0941), (5.8039, 17.5457), (6.0392, 14.8099), (6.2745, 11.9481),
        (6.3922, 9.8597), (6.7059, 8.1363), (6.9608, 7.5175), (7.3333, 6.9457),
        (7.8431, 6.274), (8.3137, 5.8627), (8.8235, 5.4168), (9.3333, 5.0617),
        (9.9216, 4.6767), (10.3725, 4.47), (10.9608, 4.2244), (11.2941, 4.13),
    ];
    validate_moody_isobar(20.00, &data, 1e-2);
}







// 
// For p0/p_ref = 30.0
// "x","y"
// 0.6667,42.3637
// 0.9412,42.3637
// 1.3529,41.8876
// 1.8235,41.4169
// 2.2941,41.4169
// 2.7059,40.9515
// 3.1765,39.5864
// 3.6471,38.7017
// 4.2353,36.1645
// 4.6275,34.5661
// 5.2549,30.8724
// 5.7843,27.5734
// 6.2353,23.5385
// 6.5882,20.3224
// 6.902,17.3486
// 7.1176,14.9782
// 7.5098,11.9481
// 7.8431,10.4329
// 8.1765,9.7489
// 8.4314,9.0074
// 8.7647,8.7072
// 9.1765,8.1363
// 9.7255,7.6893
// 10.1961,7.2669
// 
#[test]
#[ignore]
fn isobar_pref_30_00() {
    let data = vec![
        (0.6667,42.3637),
        (0.9412,42.3637),
        (1.3529,41.8876),
        (1.8235,41.4169),
        (2.2941,41.4169),
        (2.7059,40.9515),
        (3.1765,39.5864),
        (3.6471,38.7017),
        (4.2353,36.1645),
        (4.6275,34.5661),
        (5.2549,30.8724),
        (5.7843,27.5734),
        (6.2353,23.5385),
        (6.5882,20.3224),
        (6.902,17.3486),
        (7.1176,14.9782),
        (7.5098,11.9481),
        (7.8431,10.4329),
        (8.1765,9.7489),
        (8.4314,9.0074),
        (8.7647,8.7072),
        (9.1765,8.1363),
        (9.7255,7.6893),
        (10.1961,7.2669),
    ];
    validate_moody_isobar(30.00, &data, 1e-2);
}

