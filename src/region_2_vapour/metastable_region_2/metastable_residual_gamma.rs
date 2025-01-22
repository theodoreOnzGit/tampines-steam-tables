use uom::si::f64::*;
use crate::region_2_vapour::{pi_2, tau_2};

use super::METASTABLE_REGION_2_COEFFS_RES;


/// Returns the region-2 residual gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let pi = pi_2(p);
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    for coefficient in METASTABLE_REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii) * (tau - 0.5).powi(ji);
    }
    sum
}
/// Returns the region-2 residual gamma_metastable_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_tau_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    let pi = pi_2(p);
    for coefficient in METASTABLE_REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii) * f64::from(ji) * (tau - 0.5).powi(ji - 1);
    }
    sum
}

/// Returns the region-2 residual gamma_metastable_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_tau_tau_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    let pi = pi_2(p);
    for coefficient in METASTABLE_REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii) * f64::from(ji * (ji - 1)) * (tau - 0.5).powi(ji - 2);
    }
    sum
}

/// Returns the region-2 residual gamma_metastable_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_pi_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    let pi = pi_2(p);
    for coefficient in METASTABLE_REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii - 1) * f64::from(ii) * (tau - 0.5).powi(ji);
    }
    sum
}

/// Returns the region-2 residual gamma_metastable_pi_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_pi_pi_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    let pi = pi_2(p);
    for coefficient in METASTABLE_REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii - 2) * f64::from(ii * (ii - 1)) * (tau - 0.5).powi(ji);
    }
    sum
}

/// Returns the region-2 residual gamma_metastable_pi_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_pi_tau_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    let pi = pi_2(p);
    for coefficient in METASTABLE_REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii - 1) * f64::from(ii * ji) * (tau - 0.5).powi(ji - 1);
    }
    sum
}
