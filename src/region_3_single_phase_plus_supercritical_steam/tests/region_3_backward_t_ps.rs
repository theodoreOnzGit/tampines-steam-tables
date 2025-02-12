use uom::si::{f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin,  thermodynamic_temperature::kelvin};

use crate::region_3_single_phase_plus_supercritical_steam::t_ps_flash::{t_ps_3a, t_ps_3b};

#[test] 
pub fn t3a_ps_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let t_verification_kelvin = 6.282_959_869e2;


    let v_test_m3_per_kg = 
        t_ps_3a(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn t3a_ps_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.6);
    let t_verification_kelvin = 6.297_158_726e2;


    let v_test_m3_per_kg = 
        t_ps_3a(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn t3a_ps_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.0);
    let t_verification_kelvin = 7.056_880_237e2;


    let v_test_m3_per_kg = 
        t_ps_3a(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}


#[test] 
pub fn t3b_ps_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let t_verification_kelvin = 6.401_176_443e2;


    let v_test_m3_per_kg = 
        t_ps_3b(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn t3b_ps_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.5);
    let t_verification_kelvin = 7.163_687_517e2;


    let v_test_m3_per_kg = 
        t_ps_3b(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn t3b_ps_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);
    let t_verification_kelvin = 8.474_332_825e2;


    let t_test_kelvin = 
        t_ps_3b(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}


