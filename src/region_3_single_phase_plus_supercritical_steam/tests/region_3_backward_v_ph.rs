use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, specific_volume::cubic_meter_per_kilogram};

use crate::region_3_single_phase_plus_supercritical_steam::{v_ph_3a, v_ph_3b};
#[test] 
pub fn v3a_ph_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let v_verification_m3_per_kg = 1.749903962e-3;


    let v_test_m3_per_kg = 
        v_ph_3a(p, h)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3a_ph_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2000.0);
    let v_verification_m3_per_kg = 1.908139035e-3;


    let v_test_m3_per_kg = 
        v_ph_3a(p, h)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3a_ph_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2100.0);
    let v_verification_m3_per_kg = 1.676229776e-3;


    let v_test_m3_per_kg = 
        v_ph_3a(p, h)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}


#[test] 
pub fn v3b_ph_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2500.0);
    let v_verification_m3_per_kg = 6.670_547_043e-3;


    let v_test_m3_per_kg = 
        v_ph_3b(p, h)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3b_ph_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2400.0);
    let v_verification_m3_per_kg = 2.801_244_590e-3;


    let v_test_m3_per_kg = 
        v_ph_3b(p, h)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
#[test] 
pub fn v3b_ph_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let v_verification_m3_per_kg = 2.404_234_998e-3;


    let v_test_m3_per_kg = 
        v_ph_3b(p, h)
        .get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        v_verification_m3_per_kg,
        v_test_m3_per_kg,
        max_relative=1e-8
        );
}
