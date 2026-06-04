use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::available_energy::kilojoule_per_kilogram;

use crate::interfaces::functional_programming::pt_flash_eqm::h_tp_eqm_two_phase;
use crate::region_4_vap_liq_equilibrium::{sat_pressure_4, sat_temp_4, tsat_hs_4, w_ps_eqm_region4_finite_diff_vol, w_ps_eqm_region4_kieffer};

#[test]
pub fn sat_pressure_test_1(){

    let ref_p_sat_mpa = 0.353658941e-2;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);

    let p_sat_test_mpa = sat_pressure_4(t)
        .get::<megapascal>();

    approx::assert_relative_eq!(
        ref_p_sat_mpa,
        p_sat_test_mpa,
        max_relative=1e-8
        );
}
#[test]
pub fn sat_pressure_test_2(){

    let ref_p_sat_mpa = 0.263889776e1;

    let t = ThermodynamicTemperature::new::<kelvin>(500.0);

    let p_sat_test_mpa = sat_pressure_4(t)
        .get::<megapascal>();

    approx::assert_relative_eq!(
        ref_p_sat_mpa,
        p_sat_test_mpa,
        max_relative=1e-8
        );
}
#[test]
pub fn sat_pressure_test_3(){

    let ref_p_sat_mpa = 0.123443146e2;

    let t = ThermodynamicTemperature::new::<kelvin>(600.0);

    let p_sat_test_mpa = sat_pressure_4(t)
        .get::<megapascal>();

    approx::assert_relative_eq!(
        ref_p_sat_mpa,
        p_sat_test_mpa,
        max_relative=1e-8
        );
}
#[test]
pub fn sat_temp_test_1(){

    let ref_t_sat_kelvin = 0.372755919e3;

    let p = Pressure::new::<megapascal>(0.1);

    let t_sat_test_kelvin = sat_temp_4(p)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t_sat_kelvin,
        t_sat_test_kelvin,
        max_relative=1e-8
        );
}
#[test]
pub fn sat_temp_test_2(){

    let ref_t_sat_kelvin = 0.453035632e3;

    let p = Pressure::new::<megapascal>(1.0);

    let t_sat_test_kelvin = sat_temp_4(p)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t_sat_kelvin,
        t_sat_test_kelvin,
        max_relative=1e-8
        );
}
#[test]
pub fn sat_temp_test_3(){

    let ref_t_sat_kelvin = 0.584149488e3;

    let p = Pressure::new::<megapascal>(10.0);

    let t_sat_test_kelvin = sat_temp_4(p)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t_sat_kelvin,
        t_sat_test_kelvin,
        max_relative=1e-8
        );
}


#[test]
pub fn hs_backward_sat_temp_test_1(){

    let ref_t_sat_kelvin = 3.468_475_498e2;

    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.3);


    let t_sat_test_kelvin = tsat_hs_4(h,s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t_sat_kelvin,
        t_sat_test_kelvin,
        max_relative=1e-8
        );
}


#[test]
pub fn hs_backward_sat_temp_test_2(){

    let ref_t_sat_kelvin = 4.251_373_305e2;

    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2400.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.0);


    let t_sat_test_kelvin = tsat_hs_4(h,s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t_sat_kelvin,
        t_sat_test_kelvin,
        max_relative=1e-8
        );
}
#[test]
pub fn hs_backward_sat_temp_test_3(){

    let ref_t_sat_kelvin = 5.225_579_013e2;

    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2500.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.5);


    let t_sat_test_kelvin = tsat_hs_4(h,s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t_sat_kelvin,
        t_sat_test_kelvin,
        max_relative=1e-8
        );
}



// test data from:
// Kieffer, S. W. (1977). Sound speed in liquid‐gas mixtures: 
// Water‐air and water‐steam. 
// Journal of Geophysical research, 82(20), 2895-2904.
//
// Fig 9.
// 1 bar saturation pressure 
//
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
// "x","y"
// 0.00001185,1.23637487
// 0.00002393,1.23637487
// 0.00007122,1.20756707
// 0.00022252,1.52862283
// 0.00054556,2.07685484
// 0.00133757,3.4074007
// 0.00276761,5.72373055
// 0.00695193,12.1709431
// 0.0143845,23.00248436
// 0.0438637,59.06511556
// 0.12139772,128.59231566
// 0.35266992,248.83097247
// 0.76596782,398.7333797


use uom::si::pressure::bar;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;

#[test]
pub fn w_px_eqm_1_bar(){

    let p = Pressure::new::<bar>(1.0);

    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001185, 1.23637487),
        (0.00002393, 1.23637487),
        (0.00007122, 1.20756707),
        (0.00022252, 1.52862283),
        (0.00054556, 2.07685484),
        (0.00133757, 3.40740070),
        (0.00276761, 5.72373055),
        (0.00695193, 12.17094310),
        (0.01438450, 23.00248436),
        (0.04386370, 59.06511556),
        (0.12139772, 128.59231566),
        (0.35266992, 248.83097247),
        (0.76596782, 398.73337970),
    ];

        for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
            // interpolate entropy at quality x
            // compute equilibrium speed of sound
            let t_sat = sat_temp_4(p);
            let h_liq = h_tp_1(t_sat, p);
            let h_vap = h_tp_2(t_sat, p);

            // we need to find the correct enthalpy
            // so we can find the entropy
            let h = *x * h_vap + (1.0 - x) * h_liq;
            let s = s_ph_eqm(p, h);
            let w_test = w_ps_eqm_region4_kieffer(p, s);
            dbg!(&(x,w_test,w_expected));
            // assert within tolerance
            approx::assert_abs_diff_eq!(
                w_test.get::<meter_per_second>().log10(),
                w_expected.log10(),
                epsilon=0.1
            );
        }
    
}
#[test]
pub fn w_px_eqm_1_bar_finite_diff_vol(){

    let p = Pressure::new::<bar>(1.0);

    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001185, 1.23637487),
        (0.00002393, 1.23637487),
        (0.00007122, 1.20756707),
        (0.00022252, 1.52862283),
        (0.00054556, 2.07685484),
        (0.00133757, 3.40740070),
        (0.00276761, 5.72373055),
        (0.00695193, 12.17094310),
        (0.01438450, 23.00248436),
        (0.04386370, 59.06511556),
        (0.12139772, 128.59231566),
        (0.35266992, 248.83097247),
        (0.76596782, 398.73337970),
    ];

        for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
            // interpolate entropy at quality x
            // compute equilibrium speed of sound
            let t_sat = sat_temp_4(p);
            let h_liq = h_tp_1(t_sat, p);
            let h_vap = h_tp_2(t_sat, p);

            // we need to find the correct enthalpy
            // so we can find the entropy
            let h = *x * h_vap + (1.0 - x) * h_liq;
            let s = s_ph_eqm(p, h);
            let w_test = w_ps_eqm_region4_finite_diff_vol(p, s);
            dbg!(&(x,w_test,w_expected));
            // assert within tolerance
            approx::assert_abs_diff_eq!(
                w_test.get::<meter_per_second>().log10(),
                w_expected.log10(),
                epsilon=0.1
            );
        }
    
}

// 5 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
#[test]
pub fn w_px_eqm_5_bar(){
    let p = Pressure::new::<bar>(5.0);
    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001306, 4.73988624),
        (0.00003049, 4.73988624),
        (0.00011018, 4.62944587),
        (0.00027676, 4.96873384),
        (0.00071225, 5.72373055),
        (0.00143845, 6.59344865),
        (0.00344225, 9.61468707),
        (0.00729723, 14.69723015),
        (0.01438450, 23.55123315),
        (0.02905079, 42.46062164),
        (0.07847600, 90.28828407),
        (0.17462454, 166.66446149),
        (0.61584821, 354.39538240),
        (0.90760052, 427.95619490),
    ];
    for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
        let t_sat = sat_temp_4(p);
        let h_liq = h_tp_1(t_sat, p);
        let h_vap = h_tp_2(t_sat, p);
        let h = *x * h_vap + (1.0 - x) * h_liq;
        let s = s_ph_eqm(p, h);
        let w_test = w_ps_eqm_region4_finite_diff_vol(p, s);
        dbg!(&(x, w_test, w_expected));
        approx::assert_abs_diff_eq!(
            w_test.get::<meter_per_second>().log10(),
            w_expected.log10(),
            epsilon=0.1
        );
    }
}

// 10 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
#[test]
pub fn w_px_eqm_10_bar(){
    let p = Pressure::new::<bar>(10.0);
    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001274, 8.74942519),
        (0.00002905, 8.95815206),
        (0.00012438, 8.54556170),
        (0.00037019, 9.17185833),
        (0.00178909, 11.07563412),
        (0.00471720, 14.69723015),
        (0.01184907, 23.55123315),
        (0.02019590, 31.99775096),
        (0.04386370, 57.68888568),
        (0.08646535, 101.58414709),
        (0.15469408, 151.66569934),
        (0.33598183, 248.83097247),
        (0.80401316, 398.73337970),
    ];
    for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
        let t_sat = sat_temp_4(p);
        let h_liq = h_tp_1(t_sat, p);
        let h_vap = h_tp_2(t_sat, p);
        let h = *x * h_vap + (1.0 - x) * h_liq;
        let s = s_ph_eqm(p, h);
        let w_test = w_ps_eqm_region4_finite_diff_vol(p, s);
        dbg!(&(x, w_test, w_expected));
        approx::assert_abs_diff_eq!(
            w_test.get::<meter_per_second>().log10(),
            w_expected.log10(),
            epsilon=0.1
        );
    }
}

// 50 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
#[test]
pub fn w_px_eqm_50_bar(){
    let p = Pressure::new::<bar>(50.0);
    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001338, 32.76109143),
        (0.00003442, 32.76109143),
        (0.00011288, 31.99775096),
        (0.00049515, 31.99775096),
        (0.00162378, 31.25219648),
        (0.00471720, 35.16212272),
        (0.01509897, 47.77282101),
        (0.03613223, 64.90627559),
        (0.07122486, 90.28828407),
        (0.12742750, 134.80091272),
        (0.23930257, 206.05975260),
        (0.47171991, 293.47883864),
        (0.84394820, 398.73337970),
    ];
    for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
        let t_sat = sat_temp_4(p);
        let h_liq = h_tp_1(t_sat, p);
        let h_vap = h_tp_2(t_sat, p);
        let h = *x * h_vap + (1.0 - x) * h_liq;
        let s = s_ph_eqm(p, h);
        let w_test = w_ps_eqm_region4_finite_diff_vol(p, s);
        dbg!(&(x, w_test, w_expected));
        approx::assert_abs_diff_eq!(
            w_test.get::<meter_per_second>().log10(),
            w_expected.log10(),
            epsilon=0.1
        );
    }
}

// 100 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
#[test]
pub fn w_px_eqm_100_bar(){
    let p = Pressure::new::<bar>(100.0);
    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001214, 66.45468401),
        (0.00006011, 68.04003136),
        (0.00030494, 68.04003136),
        (0.00127427, 68.04003136),
        (0.00765968, 68.04003136),
        (0.02171911, 78.37868128),
        (0.06464372, 101.58414709),
        (0.13055379, 138.01673225),
        (0.23930257, 206.05975260),
        (0.47171991, 307.64836227),
        (0.82373871, 408.24558968),
    ];
    for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
        let t_sat = sat_temp_4(p);
        let h_liq = h_tp_eqm_two_phase(t_sat, p, 0.0);
        let h_vap = h_tp_eqm_two_phase(t_sat, p, 1.0);
        let h = *x * h_vap + (1.0 - x) * h_liq;
        let s = s_ph_eqm(p, h);
        let w_test = w_ps_eqm_region4_finite_diff_vol(p, s);
        dbg!(&(x, w_test, w_expected));
        approx::assert_abs_diff_eq!(
            w_test.get::<meter_per_second>().log10(),
            w_expected.log10(),
            epsilon=0.2
        );
    }
}

// 200 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
#[test]
pub fn w_px_eqm_200_bar(){
    let p = Pressure::new::<bar>(200.0);
    let quality_vs_speed_of_sound_meter_per_s: Vec<(f64, f64)> = vec![
        (0.00001157, 144.68035187),
        (0.00004952, 144.68035187),
        (0.00028355, 144.68035187),
        (0.00127427, 148.13185595),
        (0.00712249, 148.13185595),
        (0.05072980, 151.66569934),
        (0.38857395, 206.05975260),
        (0.80401316, 254.76710068),
    ];
    for (x, w_expected) in quality_vs_speed_of_sound_meter_per_s.iter() {
        let t_sat = sat_temp_4(p);
        let h_liq = h_tp_eqm_two_phase(t_sat, p, 0.0);
        let h_vap = h_tp_eqm_two_phase(t_sat, p, 1.0);
        let h = *x * h_vap + (1.0 - x) * h_liq;
        let s = s_ph_eqm(p, h);
        let w_test = w_ps_eqm_region4_finite_diff_vol(p, s);
        dbg!(&(x, w_test, w_expected));
        approx::assert_abs_diff_eq!(
            w_test.get::<meter_per_second>().log10(),
            w_expected.log10(),
            epsilon=0.2
        );
    }
}
