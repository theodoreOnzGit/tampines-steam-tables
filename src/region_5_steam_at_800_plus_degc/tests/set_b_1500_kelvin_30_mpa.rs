use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::temperature_coefficient::per_kelvin;
use uom::si::velocity::meter_per_second;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;

use crate::region_5_steam_at_800_plus_degc::*;

#[test] 
pub fn specific_vol_regression_set_b(){
    let ref_vol_m3_per_kg = 0.230761299e-1;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let specific_vol_test_m3_per_kg = 
        v_tp_5(t, p).get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        ref_vol_m3_per_kg,
        specific_vol_test_m3_per_kg,
        max_relative=1e-8);

    
}

#[test] 
pub fn specific_enthalpy_regression_set_b(){
    let ref_enthalpy_kj_per_kg = 0.516723514e4;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let specific_enthalpy_test_kj_per_kg = 
        h_tp_5(t, p).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_enthalpy_kj_per_kg,
        specific_enthalpy_test_kj_per_kg,
        max_relative=1e-9);

    
}

#[test] 
pub fn specific_internal_energy_regression_set_b(){
    let ref_internal_energy_kj_per_kg = 0.447495124e4;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let specific_internal_energy_test_kj_per_kg = 
        u_tp_5(t, p).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_internal_energy_kj_per_kg,
        specific_internal_energy_test_kj_per_kg,
        max_relative=1e-9);

    
}


#[test] 
pub fn specific_entropy_regression_set_b(){
    let ref_entropy_kj_per_kg_kelvin = 0.772970133e1;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let specific_entropy_test_kj_per_kg_kelvin = 
        s_tp_5(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_entropy_kj_per_kg_kelvin,
        specific_entropy_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn cp_regression_set_b(){
    let ref_cp_kj_per_kg_kelvin = 0.272724317e1;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let cp_test_kj_per_kg_kelvin = 
        cp_tp_5(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cp_kj_per_kg_kelvin,
        cp_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn cv_regression_set_b(){
    let ref_cv_kj_per_kg_kelvin = 0.219274829e1;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let cv_test_kj_per_kg_kelvin = 
        cv_tp_5(t, p).get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cv_kj_per_kg_kelvin,
        cv_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn speed_of_sound_regression_set_b(){
    let ref_speed_of_sound_kj_per_kg_kelvin = 0.928548002e3;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let specific_speed_of_sound_test_kj_per_kg_kelvin = 
        w_tp_5(t, p).get::<meter_per_second>();

    approx::assert_relative_eq!(
        ref_speed_of_sound_kj_per_kg_kelvin,
        specific_speed_of_sound_test_kj_per_kg_kelvin,
        max_relative=1e-8);

    
}


#[test] 
pub fn isentropic_exponent_regression_set_b(){
    let ref_kappa = 0.124544481e1;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let tested_kappa = 
        kappa_tp_5(t, p).get::<ratio>();

    approx::assert_relative_eq!(
        ref_kappa,
        tested_kappa,
        max_relative=1e-8);

    
}


#[test] 
pub fn isobaric_cubic_expansion_coeff_regression_set_b(){
    let ref_alpha_v = 0.716950754e-3;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let tested_alpha_v = 
        alpha_v_tp_5(t, p).get::<per_kelvin>();

    approx::assert_relative_eq!(
        ref_alpha_v,
        tested_alpha_v,
        max_relative=1e-8);

    
}


#[test] 
pub fn isothermal_compressibility_coeff_regression_set_b(){
    let ref_kappa_t = 0.332881253e-1;
    let t = ThermodynamicTemperature::new::<kelvin>(1500.0);
    let p = Pressure::new::<megapascal>(30.0);

    let tested_kappa_t_inverse = 
        kappa_t_tp_5(t, p).recip().get::<megapascal>();
    let tested_kappa_t = tested_kappa_t_inverse.recip();

    approx::assert_relative_eq!(
        ref_kappa_t,
        tested_kappa_t,
        max_relative=1e-8);

    
}

