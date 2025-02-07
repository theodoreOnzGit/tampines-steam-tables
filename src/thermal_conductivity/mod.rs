use std::ops::Index;

use uom::si::{f64::*, ratio::ratio};

use crate::constants::{p_crit_water, rho_crit_water, t_crit_water};

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

const LAMBDA_1_COEFFS_NI1: [[f64; 2]; 5] = [
    [1.0,  0.160_397_357e1],
    [2.0,  0.233_771_842e1],
    [3.0,  0.219_650_529e1],
    [4.0,  -0.121_051_378e1],
    [5.0,  -0.272_033_700e1],
];

const LAMBDA_1_COEFFS_NI2: [[f64; 2]; 5] = [
    [1.0,  -0.646_013_523],
    [2.0,  -0.278_843_778e1],
    [3.0,  -0.454_580_785e1],
    [4.0,  0.160_812_989e1],
    [5.0,  0.457_586_331e1],
];

const LAMBDA_1_COEFFS_NI3: [[f64; 2]; 5] = [
    [1.0,  0.111_443_906],
    [2.0,  0.153_616_167e1],
    [3.0,  0.355_777_244e1],
    [4.0,  -0.621_178_141],
    [5.0,  -0.318_369_245e1],
];
const LAMBDA_1_COEFFS_NI4: [[f64; 2]; 5] = [
    [1.0,  0.102_997_357],
    [2.0,  -0.463_045_512],
    [3.0,  -0.140_944_978e1],
    [4.0,  0.716_373_224e-1],
    [5.0,  0.111_683_480e1],
];
const LAMBDA_1_COEFFS_NI5: [[f64; 2]; 5] = [
    [1.0,  -0.504_123_634e-1],
    [2.0,  0.832_827_019e-1],
    [3.0,  0.275_418_278],
    [4.0,  0.0],
    [5.0,  -0.192_683_050],
];
const LAMBDA_1_COEFFS_NI6: [[f64; 2]; 5] = [
    [1.0,  0.609_859_258e-2],
    [2.0,  -0.719_201_245e-2],
    [3.0,  -0.205_938_816e-1],
    [4.0,  0.0],
    [5.0,  0.129_138_420e-1],
];
pub(crate) fn lambda_1(t: ThermodynamicTemperature,
    rho: MassDensity) -> f64 {
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();
    let rho_c = rho_crit_water();
    let delta_f64: f64 = (rho/rho_c).get::<ratio>();

    fn inner_sum_over_all_j(i: usize, delta: f64) -> f64{

        let ni1 = LAMBDA_1_COEFFS_NI1.index(i-1)[1];
        let ni2 = LAMBDA_1_COEFFS_NI2.index(i-1)[1];
        let ni3 = LAMBDA_1_COEFFS_NI3.index(i-1)[1];
        let ni4 = LAMBDA_1_COEFFS_NI4.index(i-1)[1];
        let ni5 = LAMBDA_1_COEFFS_NI5.index(i-1)[1];
        let ni6 = LAMBDA_1_COEFFS_NI6.index(i-1)[1];

        let inner_sum: f64 = 
            ni1 * (delta - 1.0).powf(1.0 -  1.0)
            + ni2 * (delta - 1.0).powf(2.0 - 1.0)
            + ni3 * (delta - 1.0).powf(3.0 - 1.0)
            + ni4 * (delta - 1.0).powf(4.0 - 1.0)
            + ni5 * (delta - 1.0).powf(5.0 - 1.0)
            + ni6 * (delta - 1.0).powf(6.0 - 1.0);

        return inner_sum;

    }

    
    let mut exponent: f64 = 0.0;

    // I know I should be using a for loop, but this makes coding 
    // easier to visualise given the nested sum.. lol
    let i = 1;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 2;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 3;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 4;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 5;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);


    exponent.exp()

}

#[cfg(test)]
mod tests;
