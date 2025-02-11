
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::thermodynamic_temperature::kelvin;

use crate::region_2_vapour::subregion_2a::t_ps_2a;
use crate::region_2_vapour::subregion_2b::t_ps_2b;
use crate::region_2_vapour::subregion_2c::t_ps_2c;

#[test]
pub fn subregion_2a_t_ps_test_1(){
    let p = Pressure::new::<megapascal>(0.1);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.5);
    let t_expected_kelvin = 0.399_517_097e3;

    let t_test = t_ps_2a(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2a_t_ps_test_2(){
    let p = Pressure::new::<megapascal>(0.1);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(8.0);
    let t_expected_kelvin = 0.514_127_081e3;

    let t_test = t_ps_2a(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}

#[test]
pub fn subregion_2a_t_ps_test_3(){
    let p = Pressure::new::<megapascal>(2.5);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(8.0);
    let t_expected_kelvin = 0.103_984_917e4;

    let t_test = t_ps_2a(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2b_t_ps_test_1(){
    let p = Pressure::new::<megapascal>(8.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.0);
    let t_expected_kelvin = 0.600_484_040e3;

    let t_test = t_ps_2b(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2b_t_ps_test_2(){
    let p = Pressure::new::<megapascal>(8.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.5);
    let t_expected_kelvin = 0.106_495_556e4;

    let t_test = t_ps_2b(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}

#[test]
pub fn subregion_2b_t_ps_test_3(){
    let p = Pressure::new::<megapascal>(90.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.0);
    let t_expected_kelvin = 0.103_801_126e4;

    let t_test = t_ps_2b(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2c_t_ps_test_1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.75);
    let t_expected_kelvin = 0.697_992_849e3;

    let t_test = t_ps_2c(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2c_t_ps_test_2(){
    let p = Pressure::new::<megapascal>(80.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.25);
    let t_expected_kelvin = 0.854_011_484e3;

    let t_test = t_ps_2c(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}

#[test]
pub fn subregion_2c_t_ps_test_3(){
    let p = Pressure::new::<megapascal>(80.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.75);
    let t_expected_kelvin = 0.949_017_998e3;

    let t_test = t_ps_2c(p, s);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}
