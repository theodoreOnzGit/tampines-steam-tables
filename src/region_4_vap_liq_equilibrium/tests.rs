use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::available_energy::kilojoule_per_kilogram;

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
                epsilon=0.3
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
                epsilon=0.3
            );
        }
    
}

// 5 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
// "x","y"
// 0.00001213,27.97276489
// 0.00002221,27.5335711
// 0.00005306,27.5335711
// 0.00011506,27.5335711
// 0.0002556,28.41896436
// 0.00056781,30.759021
// 0.00114505,33.82280502
// 0.00266975,41.54862932
// 0.00443669,50.23783646
// 0.00737305,59.79051613
// 0.0122528,74.61948265
// 0.01761158,93.12626067
// 0.02721911,116.22300395
// 0.04415277,152.10047986
// 0.06346306,183.90977412
// 0.10546522,233.18338074
// 0.17107765,310.03372243
// 0.30569907,393.09880013
// 0.49588178,467.84618539
// 0.80438171,548.0644279

// 10 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
// "x","y"
// 0.00001184,42.88470288
// 0.00002017,42.21138026
// 0.00003782,42.21138026
// 0.00007093,42.21138026
// 0.00013962,41.54862932
// 0.00028846,42.88470288
// 0.000541,44.26374038
// 0.00103946,46.41588834
// 0.00156812,48.67267591
// 0.00254368,51.85332767
// 0.00433067,58.85175931
// 0.006857,68.94264997
// 0.01139522,79.49570087
// 0.01677992,94.61173728
// 0.02354229,107.38109667
// 0.03224062,127.79951616
// 0.05229828,162.03990994
// 0.08082822,208.73132447
// 0.11069236,244.52099323
// 0.14796809,286.44726077
// 0.22868817,346.35295744
// 0.3131831,386.92684842
// 0.45015424,439.14878331
// 0.63156807,506.36928416
// 0.93000952,583.87922655

// 50 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
//
// "x","y"
// 0.00001369,104.03563698
// 0.00002568,104.03563698
// 0.00004702,104.03563698
// 0.00009952,102.40219794
// 0.0002056,104.03563698
// 0.00043517,102.40219794
// 0.00089905,104.03563698
// 0.00168613,102.40219794
// 0.0032397,104.03563698
// 0.00719686,114.39821393
// 0.01677992,134.01325826
// 0.03466698,159.49575946
// 0.05902125,189.82373547
// 0.08691118,222.37145504
// 0.13432327,264.65518833
// 0.21788899,340.9149512
// 0.37096055,405.73962426
// 0.54625441,475.3088985
// 0.84424906,548.0644279


// 100 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
// "x","y"
// 0.00001243,167.25060506
// 0.00002332,167.25060506
// 0.00004589,164.62464269
// 0.00009034,164.62464269
// 0.0002056,167.25060506
// 0.00052807,164.62464269
// 0.00117308,167.25060506
// 0.00242356,167.25060506
// 0.00500703,167.25060506
// 0.01167419,169.91845472
// 0.02243057,186.84335769
// 0.04206778,199.05315803
// 0.07701133,229.52222338
// 0.11902282,256.40985194
// 0.16299898,295.65850599
// 0.21788899,346.35295744
// 0.32085035,393.09880013
// 0.49588178,460.50064258
// 0.76639697,539.45939735



// 200 bar 
// x = mass fraction of steam (quality)
// y = speed of sound (m/s)
// 
// "x","y"
// 0.00001437,281.94981977
// 0.00002389,277.52299203
// 0.00006135,281.94981977
// 0.00016943,281.94981977
// 0.00044582,277.52299203
// 0.00135629,277.52299203
// 0.00340026,277.52299203
// 0.00702488,277.52299203
// 0.01451325,277.52299203
// 0.03466698,281.94981977
// 0.07517102,291.01644139
// 0.13432327,310.03372243
// 0.30569907,335.56232582
// 0.60174399,374.87213662
// 0.82407433,399.36920194

