use uom::si::{f64::*, ratio::ratio};

use crate::constants::t_crit_water;

const LAMBDA_0_COEFFS: [[f64; 2]; 5] = [
    [1.0,  0.244_322_1e-2],
    [2.0,  0.132_309_5e-1],
    [3.0,  0.677_035_7e-2],
    [4.0,  -0.345_458_6e-2],
    [5.0,  0.409_626_6e-3],
];

pub(crate) fn lambda_0(t: ThermodynamicTemperature) -> f64 {
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();

    let num = theta_f64.sqrt();

    let mut den = 0.0;

    for coeffs in LAMBDA_0_COEFFS {
        let i = coeffs[0];
        let ni = coeffs[1];

        den += ni * theta_f64.powf(1.0 - i);
    };

    return num/den;
}

#[cfg(test)]
mod tests;
