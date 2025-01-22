use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::velocity::meter_per_second;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;

use crate::region_1_subcooled_liquid::{cp_tp_1, cv_tp_1, h_tp_1, s_tp_1, u_tp_1, v_tp_1, w_tp_1};

#[test] 
pub fn specific_vol_regression_set_a(){
    let ref_vol_m3_per_kg = 0.100215168e-2;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_vol_test_m3_per_kg = 
        v_tp_1(t, p).get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        ref_vol_m3_per_kg,
        specific_vol_test_m3_per_kg,
        max_relative=1e-9);

    
}

#[test] 
pub fn specific_enthalpy_regression_set_a(){
    let ref_enthalpy_kj_per_kg = 0.115331273e3;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_enthalpy_test_kj_per_kg = 
        h_tp_1(t, p).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_enthalpy_kj_per_kg,
        specific_enthalpy_test_kj_per_kg,
        max_relative=1e-9);

    
}

#[test] 
pub fn specific_internal_energy_regression_set_a(){
    let ref_internal_energy_kj_per_kg = 0.112324818e3;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_internal_energy_test_kj_per_kg = 
        u_tp_1(t, p).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_internal_energy_kj_per_kg,
        specific_internal_energy_test_kj_per_kg,
        max_relative=1e-9);

    
}


#[test] 
pub fn specific_entropy_regression_set_a(){
    let ref_entropy_kj_per_kg_kelvin = 0.392294792;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_entropy_test_kj_per_kg_kelvin = 
        s_tp_1(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_entropy_kj_per_kg_kelvin,
        specific_entropy_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn specific_cp_regression_set_a(){
    let ref_cp_kj_per_kg_kelvin = 0.417301218e1;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_cp_test_kj_per_kg_kelvin = 
        cp_tp_1(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cp_kj_per_kg_kelvin,
        specific_cp_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn specific_cv_regression_set_a(){
    let ref_cv_kj_per_kg_kelvin = 0.412120160e1;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_cv_test_kj_per_kg_kelvin = 
        cv_tp_1(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cv_kj_per_kg_kelvin,
        specific_cv_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn specific_speed_of_sound_regression_set_a(){
    let ref_speed_of_sound_kj_per_kg_kelvin = 0.150773921e4;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_speed_of_sound_test_kj_per_kg_kelvin = 
        w_tp_1(t, p).get::<meter_per_second>();

    approx::assert_relative_eq!(
        ref_speed_of_sound_kj_per_kg_kelvin,
        specific_speed_of_sound_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}
