use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::quantities::SpecificHeatCapacity;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;

use crate::region_1_subcooled_liquid::backward_eqn_hs_1::p_hs_1;

#[test]
pub fn p_hs_flash_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(0.001);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.0);
    let p_ref = Pressure::new::<megapascal>(9.800_980_612e-4);

    let p_test = p_hs_1(h, s);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
}

#[test]
pub fn p_hs_flash_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(90.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.0);
    let p_ref = Pressure::new::<megapascal>(9.192_954_727e1);

    let p_test = p_hs_1(h, s);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
}

#[test]
pub fn p_hs_flash_test_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1500.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.4);
    let p_ref = Pressure::new::<megapascal>(5.868_294_423e1);

    let p_test = p_hs_1(h, s);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
}
