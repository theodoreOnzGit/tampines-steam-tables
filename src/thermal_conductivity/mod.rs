use uom::si::{f64::*, ratio::ratio};

use crate::constants::{p_crit_water, t_crit_water};

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
    p: Pressure) -> f64 {
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();
    let p_c = p_crit_water();
    let delta_f64: f64 = (p/p_c).get::<ratio>();

    fn inner_sum(i: f64, delta: f64, j: f64) -> f64{

        let coeff_arr: [[f64;2]; 5];

        if j == 1.0 {
            coeff_arr = LAMBDA_1_COEFFS_NI1;
        } else if j == 2.0 {
            coeff_arr = LAMBDA_1_COEFFS_NI2;
        } else if j == 3.0 {
            coeff_arr = LAMBDA_1_COEFFS_NI3;
        } else if j == 4.0 {
            coeff_arr = LAMBDA_1_COEFFS_NI4;
        } else if j == 5.0 {
            coeff_arr = LAMBDA_1_COEFFS_NI5;
        } else {
            coeff_arr = LAMBDA_1_COEFFS_NI6;
        };

        let mut sum = 0.0;

        for coeffs in coeff_arr {
            let nij = coeffs[1];

            sum += nij * (delta - 1.0).powf(j-1.0);
        }

        sum

    }

    todo!()

}

#[cfg(test)]
mod tests;
