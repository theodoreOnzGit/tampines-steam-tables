use uom::si::{f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin, specific_volume::cubic_meter_per_kilogram};

use crate::region_3_single_phase_plus_supercritical_steam::{v_ps_3, v_ps_flash::{v_ps_3a, v_ps_3b}};
#[test] 
pub fn v3a_ps_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let v_verification_m3_per_kg = 1.733_791_463e-3;


    let v_test_m3_per_kg = 
        v_ps_3a(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3a_ps_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.6);
    let v_verification_m3_per_kg = 1.469_680_170e-3;


    let v_test_m3_per_kg = 
        v_ps_3a(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3a_ps_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.0);
    let v_verification_m3_per_kg = 1.555_893_131e-3;


    let v_test_m3_per_kg = 
        v_ps_3a(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}


#[test] 
pub fn v3b_ps_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let v_verification_m3_per_kg = 6.262_101_987e-3;


    let v_test_m3_per_kg = 
        v_ps_3b(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3b_ps_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.5);
    let v_verification_m3_per_kg = 2.332_634_294e-3;


    let v_test_m3_per_kg = 
        v_ps_3b(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3b_ps_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let v_verification_m3_per_kg = 2.449_610_757e-3;


    let v_test_m3_per_kg = 
        v_ps_3b(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}


#[test] 
pub fn v3_region_a_ps_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let v_verification_m3_per_kg = 1.733_791_463e-3;


    let v_test_m3_per_kg = 
        v_ps_3(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3_region_a_ps_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.6);
    let v_verification_m3_per_kg = 1.469_680_170e-3;


    let v_test_m3_per_kg = 
        v_ps_3(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3_region_a_ps_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.0);
    let v_verification_m3_per_kg = 1.555_893_131e-3;


    let v_test_m3_per_kg = 
        v_ps_3(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}


#[test] 
pub fn v3_reigon_b_ps_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let v_verification_m3_per_kg = 6.262_101_987e-3;


    let v_test_m3_per_kg = 
        v_ps_3(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3_reigon_b_ps_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.5);
    let v_verification_m3_per_kg = 2.332_634_294e-3;


    let v_test_m3_per_kg = 
        v_ps_3(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3_reigon_b_ps_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let v_verification_m3_per_kg = 2.449_610_757e-3;


    let v_test_m3_per_kg = 
        v_ps_3(p, s)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
