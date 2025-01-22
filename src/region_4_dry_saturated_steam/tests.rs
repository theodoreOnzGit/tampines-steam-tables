use uom::si::{f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

use crate::region_4_dry_saturated_steam::{sat_pressure_4, sat_temp_4};

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
