use super::{delta_3, tau_3, REGION_3_COEFFS};


use uom::si::f64::*;
/// Returns the region-3 phi 
/// remember, phi is dimensionless_helmholtz_free_energy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn phi_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let mut sum: f64 = 0.0;
    let tau: f64 = tau_3(t);
    let delta: f64 = delta_3(rho);
    for coefficient in REGION_3_COEFFS.iter().skip(1) {
        let ii: i32 = coefficient[0] as i32;
        let ji: i32 = coefficient[1] as i32;
        let ni: f64 = coefficient[2];
        sum += ni * delta.powi(ii) * tau.powi(ji);
    }
    sum + REGION_3_COEFFS[0][2] * delta_3(rho).ln()
}
