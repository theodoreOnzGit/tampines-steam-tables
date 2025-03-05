use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::quantities::SpecificHeatCapacity;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;

use crate::region_2_vapour::{p_hs_2a, p_hs_2b, p_hs_2c};

#[test]
pub fn subregion_2a_p_hs_test_1(){
    let p = Pressure::new::<megapascal>(1.371_012_767);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.5);

    let p_test = p_hs_2a(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2a_p_hs_test_2(){
    let p = Pressure::new::<megapascal>(1.879_743_844e-3);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.5);

    let p_test = p_hs_2a(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2a_p_hs_test_3(){
    let p = Pressure::new::<megapascal>(1.024_788_997e-1);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(4100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.5);

    let p_test = p_hs_2a(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2b_p_hs_test_1(){
    let p = Pressure::new::<megapascal>(4.793_911_442);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.0);

    let p_test = p_hs_2b(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2b_p_hs_test_2(){
    let p = Pressure::new::<megapascal>(8.395_519_209e1);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3600.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.0);

    let p_test = p_hs_2b(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2b_p_hs_test_3(){
    let p = Pressure::new::<megapascal>(7.527_161_441);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3600.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.0);

    let p_test = p_hs_2b(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2c_p_hs_test_1(){
    let p = Pressure::new::<megapascal>(9.439_202_060e1);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.1);

    let p_test = p_hs_2c(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2c_p_hs_test_2(){
    let p = Pressure::new::<megapascal>(8.414_574_124);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.8);

    let p_test = p_hs_2c(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2c_p_hs_test_3(){
    let p = Pressure::new::<megapascal>(8.376_903_879e1);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3400.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.8);

    let p_test = p_hs_2c(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
