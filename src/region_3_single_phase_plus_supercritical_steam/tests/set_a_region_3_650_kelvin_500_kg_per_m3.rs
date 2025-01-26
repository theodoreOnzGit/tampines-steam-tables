use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::temperature_coefficient::per_kelvin;
use uom::si::velocity::meter_per_second;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::region_3_single_phase_plus_supercritical_steam::{alpha_p_rho_t_3, alpha_v_rho_t_3, beta_p_rho_t_3, cp_rho_t_3, cv_rho_t_3, h_rho_t_3, kappa_rho_t_3, kappa_t_rho_t_3, p_rho_t_3, s_rho_t_3, u_rho_t_3, w_rho_t_3};

#[test] 
pub fn pressure_regression_set_a(){

    let ref_pressure_mpa = 0.255837018e2;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let pressure_test_mpa = 
        p_rho_t_3(rho, t).get::<megapascal>();

    approx::assert_relative_eq!(
        ref_pressure_mpa,
        pressure_test_mpa,
        max_relative=1e-9);
    
}

#[test] 
pub fn specific_enthalpy_regression_set_a(){
    let ref_enthalpy_kj_per_kg = 0.186343019e4;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let specific_enthalpy_test_kj_per_kg = 
      h_rho_t_3(rho, t).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_enthalpy_kj_per_kg,
        specific_enthalpy_test_kj_per_kg,
        max_relative=1e-9);

}

#[test] 
pub fn specific_internal_energy_regression_set_a(){
    let ref_internal_energy_kj_per_kg = 0.181226279e4;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let specific_internal_energy_test_kj_per_kg = 
        u_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_internal_energy_kj_per_kg,
        specific_internal_energy_test_kj_per_kg,
        max_relative=1e-8);

    
}


#[test] 
pub fn specific_entropy_regression_set_a(){
    let ref_entropy_kj_per_kg_kelvin = 0.405427273e1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let specific_entropy_test_kj_per_kg_kelvin = 
        s_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_entropy_kj_per_kg_kelvin,
        specific_entropy_test_kj_per_kg_kelvin,
        max_relative=1e-8);
    
}


#[test] 
pub fn cp_regression_set_a(){
    let ref_cp_kj_per_kg_kelvin = 0.138935717e2;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let cp_test_kj_per_kg_kelvin = 
        cp_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cp_kj_per_kg_kelvin,
        cp_test_kj_per_kg_kelvin,
        max_relative=1e-8);

}


#[test] 
pub fn cv_regression_set_a(){
    let ref_cv_kj_per_kg_kelvin = 0.319131787e1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let cv_test_kj_per_kg_kelvin = 
        cv_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cv_kj_per_kg_kelvin,
        cv_test_kj_per_kg_kelvin,
        max_relative=1e-8);
    
}


#[test] 
pub fn speed_of_sound_regression_set_a(){
    let ref_speed_of_sound_kj_per_kg_kelvin = 0.502005554e3;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let specific_speed_of_sound_test_kj_per_kg_kelvin = 
        w_rho_t_3(rho, t)
        .get::<meter_per_second>();

    approx::assert_relative_eq!(
        ref_speed_of_sound_kj_per_kg_kelvin,
        specific_speed_of_sound_test_kj_per_kg_kelvin,
        max_relative=1e-8);

}


#[test] 
pub fn isentropic_exponent_regression_set_a(){
    let ref_kappa = 0.492519765e1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let tested_kappa = 
        kappa_rho_t_3(rho, t);

    approx::assert_relative_eq!(
        ref_kappa,
        tested_kappa,
        max_relative=1e-8);

}


#[test] 
pub fn isobaric_cubic_expansion_coeff_regression_set_a(){
    let ref_alpha_v = 0.168653107e-1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let tested_alpha_v = 
        alpha_v_rho_t_3(rho, t)
        .get::<per_kelvin>();

    approx::assert_relative_eq!(
        ref_alpha_v,
        tested_alpha_v,
        max_relative=1e-8);

    
}


#[test] 
pub fn isothermal_compressibility_coeff_regression_set_a(){
    let ref_kappa_t = 0.345506956e-1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let tested_kappa_t_inverse = 
        kappa_t_rho_t_3(rho, t)
        .recip().get::<megapascal>();
    let tested_kappa_t = tested_kappa_t_inverse.recip();

    approx::assert_relative_eq!(
        ref_kappa_t,
        tested_kappa_t,
        max_relative=1e-8);

}

#[test] 
pub fn relative_pressure_coeff_regression_set_a(){
    let ref_alpha_p_per_kelvin = 0.190798153e-1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);

    let tested_alpha_p_per_kelvin = 
        alpha_p_rho_t_3(rho, t)
        .get::<per_kelvin>();

    approx::assert_relative_eq!(
        ref_alpha_p_per_kelvin,
        tested_alpha_p_per_kelvin,
        max_relative=1e-8);
    
}
#[test] 
pub fn isothermal_stress_coeff_regression_set_a(){

    let ref_beta_p_kg_per_m3 = 0.565652647e3;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(500.0);
    
    let tested_beta_p_kg_per_m3 = 
        beta_p_rho_t_3(rho, t)
        .get::<kilogram_per_cubic_meter>();
    
    approx::assert_relative_eq!(
        ref_beta_p_kg_per_m3,
        tested_beta_p_kg_per_m3,
        max_relative=1e-8);
}
