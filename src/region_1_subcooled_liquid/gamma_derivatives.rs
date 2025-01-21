
// ================    Region 1 ===================

use uom::si::f64::*;

use super::{pi_1, tau_1, REGION_1_COEFFS};

/// Returns the region-1 gamma_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_pi_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_COEFFS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += -ni * f64::from(ii) * (7.1 - pi).powi(ii - 1) * (tau - 1.222).powi(ji);
    }
    sum
}

/// Returns the region-1 gamma_pi_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_pi_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_COEFFS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * f64::from(ii * (ii - 1)) * (7.1 - pi).powi(ii - 2) * (tau - 1.222).powi(ji);
    }
    sum
}

/// Returns the region-1 gamma_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_tau_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_COEFFS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * f64::from(ji) * (7.1 - pi).powi(ii) * (tau - 1.222).powi(ji - 1);
    }
    sum
}

/// Returns the region-1 gamma_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_tau_tau_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_COEFFS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * f64::from(ji * (ji - 1)) * (7.1 - pi).powi(ii) * (tau - 1.222).powi(ji - 2);
    }
    sum
}

/// Returns the region-1 gamma_pi_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_tau_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_COEFFS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += -ni * f64::from(ii * ji) * (7.1 - pi).powi(ii - 1) * (tau - 1.222).powi(ji - 1);
    }
    sum
}
