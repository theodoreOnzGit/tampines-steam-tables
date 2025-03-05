use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::quantities::SpecificHeatCapacity;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;

use crate::region_2_vapour::p_hs_2a;

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
