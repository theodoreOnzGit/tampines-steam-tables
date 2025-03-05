use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::available_energy::kilojoule_per_kilogram;

use crate::region_4_vap_liq_equilibrium::{sat_pressure_4, sat_temp_4, tsat_hs_4};

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
