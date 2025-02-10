use super::{pi_5, tau_5, REGION_5_COEFFS_IDEAL};
use uom::si::f64::*;

/// Returns the region-5 ideal gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_5_ideal(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let pi: f64 = pi_5(p);
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    for coefficient in REGION_5_COEFFS_IDEAL {
        let ji: i32 = coefficient[0] as i32;
        let ni: f64 = coefficient[1];
        sum += ni * tau.powi(ji);
    }
    pi.ln() + sum
}

/// Returns the region-5 ideal gamma_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_tau_tau_5_ideal(t: ThermodynamicTemperature, _p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    for coefficient in REGION_5_COEFFS_IDEAL {
        let ji: f64 = coefficient[0];
        let ni: f64 = coefficient[1];
        sum += ni * ji * (ji - 1.0) * tau.powi(ji as i32 - 2);
    }
    sum
}

/// Returns the region-5 ideal gamma_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_pi_5_ideal(_t: ThermodynamicTemperature, p: Pressure) -> f64 {
    1.0 / pi_5(p)
}
/// Returns the region-5 ideal gamma_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_tau_5_ideal(t: ThermodynamicTemperature, _p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_5(t);
    for coefficient in REGION_5_COEFFS_IDEAL {
        let ji: f64 = coefficient[0];
        let ni: f64 = coefficient[1];
        sum += ni * ji * tau.powi(ji as i32 - 1);
    }
    sum
}
