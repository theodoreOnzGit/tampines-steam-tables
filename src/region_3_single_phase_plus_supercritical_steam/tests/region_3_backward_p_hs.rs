use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

use crate::region_3_single_phase_plus_supercritical_steam::p_hs_3;


#[test]
pub fn p_hs_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let p = Pressure::new::<megapascal>(2.555_703_246e1);

    let p_test = p_hs_3(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
        
}


#[test]
pub fn p_hs_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2000.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.2);
    let p = Pressure::new::<megapascal>(4.540_873_468e1);

    let p_test = p_hs_3(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
        
}


#[test]
pub fn p_hs_test_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.3);
    let p = Pressure::new::<megapascal>(6.078_123_340e1);

    let p_test = p_hs_3(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
        
}
