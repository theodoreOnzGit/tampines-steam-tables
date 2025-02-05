use uom::si::{f64::*, ratio::ratio};

use crate::constants::t_crit_water;


pub const PSI_0_COEFFS: [[f64; 2]; 4] = [
    [1.0,  0.167_752e-1],
    [2.0,  0.220_462e-1],
    [3.0,  0.636_656_4e-2],
    [4.0,  -0.241_605e-2],
];

pub(crate) fn psi_0_viscosity(t: ThermodynamicTemperature) -> f64 {
    let t_crit_water = t_crit_water();
    let theta = t/t_crit_water;
    let theta_f64 = theta.get::<ratio>();

    let mut den = 0.0;

    for coeffs in PSI_0_COEFFS {
        let i = coeffs[0];
        let ni = coeffs[1];

        den += ni * theta_f64.powf(1.0 - i);

    }

    return theta_f64.sqrt()*den.recip();


}

#[cfg(test)]
mod tests;
