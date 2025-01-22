use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::temperature_coefficient::per_kelvin;
use uom::si::velocity::meter_per_second;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;

use crate::region_2_vapour::*;

#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn specific_vol_regression_set_b(){
    let ref_vol_m3_per_kg = 0.186212297;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let specific_vol_test_m3_per_kg = 
        v_tp_2_metastable(t, p).get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        ref_vol_m3_per_kg,
        specific_vol_test_m3_per_kg,
        max_relative=1e-8);

    
}

#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn specific_enthalpy_regression_set_b(){
    let ref_enthalpy_kj_per_kg = 0.274015123e4;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let specific_enthalpy_test_kj_per_kg = 
        h_tp_2_metastable(t, p).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_enthalpy_kj_per_kg,
        specific_enthalpy_test_kj_per_kg,
        max_relative=1e-9);

    
}

#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn specific_internal_energy_regression_set_b(){
    let ref_internal_energy_kj_per_kg = 0.255393894e4;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let specific_internal_energy_test_kj_per_kg = 
        u_tp_2_metastable(t, p).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_internal_energy_kj_per_kg,
        specific_internal_energy_test_kj_per_kg,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn specific_entropy_regression_set_b(){
    let ref_entropy_kj_per_kg_kelvin = 0.650218759e1;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let specific_entropy_test_kj_per_kg_kelvin = 
        s_tp_2_metastable(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_entropy_kj_per_kg_kelvin,
        specific_entropy_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn cp_regression_set_b(){
    let ref_cp_kj_per_kg_kelvin = 0.298166443e1;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let cp_test_kj_per_kg_kelvin = 
        cp_tp_2_metastable(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cp_kj_per_kg_kelvin,
        cp_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn cv_regression_set_b(){
    let ref_cv_kj_per_kg_kelvin = 0.208622142e1;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let cv_test_kj_per_kg_kelvin = 
        cv_tp_2_metastable(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cv_kj_per_kg_kelvin,
        cv_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn speed_of_sound_regression_set_b(){
    let ref_speed_of_sound_kj_per_kg_kelvin = 0.489363295e3;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let specific_speed_of_sound_test_kj_per_kg_kelvin = 
        w_tp_2_metastable(t, p).get::<meter_per_second>();

    approx::assert_relative_eq!(
        ref_speed_of_sound_kj_per_kg_kelvin,
        specific_speed_of_sound_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn isentropic_exponent_regression_set_b(){
    let ref_kappa = 0.128603985e1;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let tested_kappa = 
        kappa_tp_2_metastable(t, p).get::<ratio>();

    approx::assert_relative_eq!(
        ref_kappa,
        tested_kappa,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn isobaric_cubic_expansion_coeff_regression_set_b(){
    let ref_alpha_v = 0.348506136e-2;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let tested_alpha_v = 
        alpha_v_tp_2_metastable(t, p).get::<per_kelvin>();

    approx::assert_relative_eq!(
        ref_alpha_v,
        tested_alpha_v,
        max_relative=1e-8);

    
}


#[test] 
//#[ignore = "metastable region 2 todo"]
pub fn isothermal_compressibility_coeff_regression_set_b(){
    let ref_kappa_t = 0.111133230e1;
    let t = ThermodynamicTemperature::new::<kelvin>(440.0);
    let p = Pressure::new::<megapascal>(1.0);

    let tested_kappa_t_inverse = 
        kappa_t_tp_2_metastable(t, p).recip().get::<megapascal>();
    let tested_kappa_t = tested_kappa_t_inverse.recip();

    approx::assert_relative_eq!(
        ref_kappa_t,
        tested_kappa_t,
        max_relative=1e-8);

    
}
