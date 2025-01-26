use super::{delta_3, tau_3, REGION_3_COEFFS};
use uom::si::f64::*;


/// Returns the region-3 phi_delta
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn phi_delta_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_3(t);
    let delta: f64 = delta_3(rho);
    for coefficient in REGION_3_COEFFS.iter().skip(1) {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * delta.powi(ii - 1) * f64::from(ii) * tau.powi(ji);
    }
    sum + REGION_3_COEFFS[0][2] / delta
}

/// Returns the region-3 phi_delta_delta
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn phi_delta_delta_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_3(t);
    let delta: f64 = delta_3(rho);
    for coefficient in REGION_3_COEFFS.iter().skip(1) {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * delta.powi(ii - 2) * f64::from(ii) * f64::from(ii - 1) * tau.powi(ji);
    }
    sum - REGION_3_COEFFS[0][2] / delta.powi(2)
}

/// Returns the region-3 phi_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn phi_tau_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_3(t);
    let delta: f64 = delta_3(rho);
    for coefficient in REGION_3_COEFFS.iter().skip(1) {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * delta.powi(ii) * f64::from(ji) * tau.powi(ji - 1);
    }
    sum
}

/// Returns the region-3 phi_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn phi_tau_tau_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_3(t);
    let delta: f64 = delta_3(rho);
    for coefficient in REGION_3_COEFFS.iter().skip(1) {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * delta.powi(ii) * f64::from(ji) * f64::from(ji - 1) * tau.powi(ji - 2);
    }
    sum
}

/// Returns the region-3 phi_delta_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn phi_delta_tau_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_3(t);
    let delta: f64 = delta_3(rho);
    for coefficient in REGION_3_COEFFS.iter().skip(1) {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * delta.powi(ii - 1) * f64::from(ii) * f64::from(ji) * tau.powi(ji - 1);
    }
    sum
}
