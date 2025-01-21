
pub const REGION_1_COEFFS: [[f64; 3]; 34] = [
    [0.0, -2.0, 0.14632971213167],
    [0.0, -1.0, -0.84548187169114],
    [0.0, 0.0, -0.37563603672040e1],
    [0.0, 1.0, 0.33855169168385e1],
    [0.0, 2.0, -0.95791963387872],
    [0.0, 3.0, 0.15772038513228],
    [0.0, 4.0, -0.16616417199501e-1],
    [0.0, 5.0, 0.81214629983568e-3],
    [1.0, -9.0, 0.28319080123804e-3],
    [1.0, -7.0, -0.60706301565874e-3],
    [1.0, -1.0, -0.18990068218419e-1],
    [1.0, 0.0, -0.32529748770505e-1],
    [1.0, 1.0, -0.21841717175414e-1],
    [1.0, 3.0, -0.52838357969930e-4],
    [2.0, -3.0, -0.47184321073267e-3],
    [2.0, 0.0, -0.30001780793026e-3],
    [2.0, 1.0, 0.47661393906987e-4],
    [2.0, 3.0, -0.44141845330846e-5],
    [2.0, 17.0, -0.72694996297594e-15],
    [3.0, -4.0, -0.31679644845054e-4],
    [3.0, 0.0, -0.28270797985312e-5],
    [3.0, 6.0, -0.85205128120103e-9],
    [4.0, -5.0, -0.22425281908000e-5],
    [4.0, -2.0, -0.65171222895601e-6],
    [4.0, 10.0, -0.14341729937924e-12],
    [5.0, -8.0, -0.40516996860117e-6],
    [8.0, -11.0, -0.12734301741641e-8],
    [8.0, -6.0, -0.17424871230634e-9],
    [21.0, -29.0, -0.68762131295531e-18],
    [23.0, -31.0, 0.14478307828521e-19],
    [29.0, -38.0, 0.26335781662795e-22],
    [30.0, -39.0, -0.11947622640071e-22],
    [31.0, -40.0, 0.18228094581404e-23],
    [32.0, -41.0, -0.93537087292458e-25],
];

// ================    Region 1 ===================

use uom::si::f64::*;

use super::{pi_1, tau_1};
/// Returns the region-1 gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
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

/// Returns the region-1 gamma_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_1(t: ThermodynamicTemperature, p: Pressure) -> f64 {
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


