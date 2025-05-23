use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::temperature_coefficient::per_kelvin;
use uom::si::velocity::meter_per_second;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::region_3_single_phase_plus_supercritical_steam::{alpha_p_rho_t_3, alpha_v_rho_t_3, beta_p_rho_t_3, cp_rho_t_3, cv_rho_t_3, h_rho_t_3, kappa_rho_t_3, kappa_t_rho_t_3, p_rho_t_3, s_rho_t_3, u_rho_t_3, w_rho_t_3};

#[test] 
pub fn pressure_regression_set_b(){

    let ref_pressure_mpa = 0.222930643e2;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let pressure_test_mpa = 
        p_rho_t_3(rho, t).get::<megapascal>();

    approx::assert_relative_eq!(
        ref_pressure_mpa,
        pressure_test_mpa,
        max_relative=1e-8);
    
}

#[test] 
pub fn specific_enthalpy_regression_set_b(){
    let ref_enthalpy_kj_per_kg = 0.237512401e4;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let specific_enthalpy_test_kj_per_kg = 
      h_rho_t_3(rho, t).get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_enthalpy_kj_per_kg,
        specific_enthalpy_test_kj_per_kg,
        max_relative=1e-8);

}

#[test] 
pub fn specific_internal_energy_regression_set_b(){
    let ref_internal_energy_kj_per_kg = 0.226365868e4;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let specific_internal_energy_test_kj_per_kg = 
        u_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        ref_internal_energy_kj_per_kg,
        specific_internal_energy_test_kj_per_kg,
        max_relative=1e-8);

    
}


#[test] 
pub fn specific_entropy_regression_set_b(){
    let ref_entropy_kj_per_kg_kelvin = 0.485438792e1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let specific_entropy_test_kj_per_kg_kelvin = 
        s_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_entropy_kj_per_kg_kelvin,
        specific_entropy_test_kj_per_kg_kelvin,
        max_relative=1e-8);
    
}


#[test] 
pub fn cp_regression_set_b(){
    let ref_cp_kj_per_kg_kelvin = 0.446579342e2;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let cp_test_kj_per_kg_kelvin = 
        cp_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cp_kj_per_kg_kelvin,
        cp_test_kj_per_kg_kelvin,
        max_relative=1e-8);

}


#[test] 
pub fn cv_regression_set_b(){
    let ref_cv_kj_per_kg_kelvin = 0.404118076e1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let cv_test_kj_per_kg_kelvin = 
        cv_rho_t_3(rho, t)
        .get::<kilojoule_per_kilogram_kelvin>();

    approx::assert_relative_eq!(
        ref_cv_kj_per_kg_kelvin,
        cv_test_kj_per_kg_kelvin,
        max_relative=1e-8);
    
}


#[test] 
pub fn speed_of_sound_regression_set_b(){
    let ref_speed_of_sound_kj_per_kg_kelvin = 0.383444594e3;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let specific_speed_of_sound_test_kj_per_kg_kelvin = 
        w_rho_t_3(rho, t)
        .get::<meter_per_second>();

    approx::assert_relative_eq!(
        ref_speed_of_sound_kj_per_kg_kelvin,
        specific_speed_of_sound_test_kj_per_kg_kelvin,
        max_relative=1e-8);

}


#[test] 
pub fn isentropic_exponent_regression_set_b(){
    let ref_kappa = 0.131906278e1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let tested_kappa = 
        kappa_rho_t_3(rho, t);

    approx::assert_relative_eq!(
        ref_kappa,
        tested_kappa,
        max_relative=1e-8);

}


#[test] 
pub fn isobaric_cubic_expansion_coeff_regression_set_b(){
    let ref_alpha_v = 0.685312229e-1;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let tested_alpha_v = 
        alpha_v_rho_t_3(rho, t)
        .get::<per_kelvin>();

    approx::assert_relative_eq!(
        ref_alpha_v,
        tested_alpha_v,
        max_relative=1e-8);

    
}


#[test] 
pub fn isothermal_compressibility_coeff_regression_set_b(){
    let ref_kappa_t = 0.375798565;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

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
pub fn relative_pressure_coeff_regression_set_b(){
    let ref_alpha_p_per_kelvin = 0.818019386e-2;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);

    let tested_alpha_p_per_kelvin = 
        alpha_p_rho_t_3(rho, t)
        .get::<per_kelvin>();

    approx::assert_relative_eq!(
        ref_alpha_p_per_kelvin,
        tested_alpha_p_per_kelvin,
        max_relative=1e-8);
    
}
#[test] 
pub fn isothermal_stress_coeff_regression_set_b(){

    let ref_beta_p_kg_per_m3 = 0.238728962e2;
    let t = ThermodynamicTemperature::new::<kelvin>(650.0);
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(200.0);
    
    let tested_beta_p_kg_per_m3 = 
        beta_p_rho_t_3(rho, t)
        .get::<kilogram_per_cubic_meter>();
    
    approx::assert_relative_eq!(
        ref_beta_p_kg_per_m3,
        tested_beta_p_kg_per_m3,
        max_relative=1e-8);
}

