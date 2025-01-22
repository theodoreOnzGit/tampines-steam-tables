use uom::si::f64::*;

use crate::region_2_vapour::{pi_2, tau_2};

use super::METASTABLE_REGION_2_COEFFS_IDEAL;


/// Returns the region-2 ideal gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_2_ideal(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let pi = pi_2(p);
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    for coefficient in METASTABLE_REGION_2_COEFFS_IDEAL {
        let ji = coefficient[0];
        let ni = coefficient[1];
        sum += ni * tau.powf(ji);
    }
    pi.ln() + sum
}


/// Returns the region-2 ideal gamma_metastable_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_tau_2_ideal(t: ThermodynamicTemperature, _p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    for coefficient in METASTABLE_REGION_2_COEFFS_IDEAL {
        let ji = coefficient[0];
        let ni = coefficient[1];
        sum += ni * ji * tau.powf(ji - 1.0);
    }
    sum
}

/// Returns the region-2 ideal gamma_metastable_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_tau_tau_2_ideal(t: ThermodynamicTemperature, _p: Pressure) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    for coefficient in METASTABLE_REGION_2_COEFFS_IDEAL {
        let ji = coefficient[0];
        let ni = coefficient[1];
        sum += ni * ji * (ji - 1.0) * tau.powf(ji - 2.0);
    }
    sum
}

/// Returns the region-2 ideal gamma_metastable_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_metastable_pi_2_ideal(_t: ThermodynamicTemperature, p: Pressure) -> f64 {
    1.0 / pi_2(p)
}
