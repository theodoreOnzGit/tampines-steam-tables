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


#[test]
pub fn p_hs_test_4(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2400.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.7);
    let p = Pressure::new::<megapascal>(6.363_924_887e1);

    let p_test = p_hs_3(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
        
}


#[test]
pub fn p_hs_test_5(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.1);
    let p = Pressure::new::<megapascal>(3.434_999_263e1);

    let p_test = p_hs_3(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
        
}


#[test]
pub fn p_hs_test_6(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let p = Pressure::new::<megapascal>(8.839_043_281e1);

    let p_test = p_hs_3(h, s);

    approx::assert_relative_eq!(
        p.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
        
}
