use super::{pi_2, tau_2, REGION_2_COEFFS_RES};
use uom::si::f64::*;


/// Returns the region-2 residual gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_2_res(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let pi = pi_2(p);
    let mut sum: f64 = 0.0;
    let tau = tau_2(t);
    for coefficient in REGION_2_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii) * (tau - 0.5).powi(ji);
    }
    sum
}
