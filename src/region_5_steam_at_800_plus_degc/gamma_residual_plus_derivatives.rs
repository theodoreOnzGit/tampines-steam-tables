use super::{pi_5, tau_5, REGION_5_COEFFS_RES};
use uom::si::f64::*;


/// Returns the region-2 residual gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_5_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let pi: f64 = pi_5(p);
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    for coefficient in REGION_5_COEFFS_RES {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * pi.powi(ii) * tau.powi(ji);
    }
    sum
}

/// Returns the region-5 residual gamma_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_tau_tau_5_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    let pi: f64 = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii: i32 = coefficient[0] as i32;
        let ji: f64 = coefficient[1];
        let ni: f64 = coefficient[2];
        sum += ni * pi.powi(ii) * ji * (ji - 1.0) * tau.powi(ji as i32 - 2);
    }
    sum
}

/// Returns the region-5 residual gamma_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_pi_5_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    let pi: f64 = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii: f64 = coefficient[0];
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * ii * pi.powi(ii as i32 - 1) * tau.powi(ji);
    }
    sum
}

/// Returns the region-5 residual gamma_pi_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_pi_pi_5_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    let pi: f64 = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii: f64 = coefficient[0];
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * ii * (ii - 1.0) * pi.powi(ii as i32 - 2) * tau.powi(ji);
    }
    sum
}

/// Returns the region-5 residual gamma_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_tau_5_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    let pi: f64 = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii: i32 = coefficient[0] as i32;
        let ji: f64 = coefficient[1];
        let ni: f64 = coefficient[2];
        sum += ni * pi.powi(ii) * ji * tau.powi(ji as i32 - 1);
    }
    sum
}
/// Returns the region-5 residual gamma_pi_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_pi_tau_5_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    let pi: f64 = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii: f64 = coefficient[0];
        let ji: f64 = coefficient[1];
        let ni: f64 = coefficient[2];
        sum += ni * ii * pi.powi(ii as i32 - 1) * ji * tau.powi(ji as i32 - 1);
    }
    sum
}
