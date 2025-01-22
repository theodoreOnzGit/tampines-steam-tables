
// ================    Region 1 ===================

use uom::si::f64::*;

use super::{pi_1, tau_1, REGION_1_COEFFS};
/// Returns the region-1 gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn gamma_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_COEFFS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * (7.1 - pi).powi(ii) * (tau - 1.222).powi(ji);
    }
    sum
}



