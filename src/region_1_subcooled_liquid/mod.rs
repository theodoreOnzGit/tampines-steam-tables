
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

pub const REGION_1_BACK_COEFFS_PH: [[f64; 3]; 20] = [
    [0.0, 0.0, -0.23872489924521e+3],
    [0.0, 1.0, 0.40421188637945e+3],
    [0.0, 2.0, 0.11349746881718e+3],
    [0.0, 6.0, -0.58457616048039e+1],
    [0.0, 22.0, -0.15285482413140e-3],
    [0.0, 32.0, -0.10866707695377e-5],
    [1.0, 0.0, -0.13391744872602e+2],
    [1.0, 1.0, 0.43211039183559e+2],
    [1.0, 2.0, -0.54010067170506e+2],
    [1.0, 3.0, 0.30535892203916e+2],
    [1.0, 4.0, -0.65964749423638e+1],
    [1.0, 10.0, 0.93965400878363e-2],
    [1.0, 32.0, 0.11573647505340e-6],
    [2.0, 10.0, -0.25858641282073e-4],
    [2.0, 32.0, -0.40644363084799e-8],
    [3.0, 10.0, 0.66456186191635e-7],
    [3.0, 32.0, 0.80670734103027e-10],
    [4.0, 32.0, -0.93477771213947e-12],
    [5.0, 32.0, 0.58265442020601e-14],
    [6.0, 32.0, -0.15020185953503e-16],
];
const REGION_1_BACK_COEFFS_PS: [[f64; 3]; 20] = [
    [0.0, 0.0, 0.17478268058307e+03],
    [0.0, 1.0, 0.34806930892873e+02],
    [0.0, 2.0, 0.65292584978455e+01],
    [0.0, 3.0, 0.33039981775489],
    [0.0, 11.0, -0.19281382923196e-06],
    [0.0, 31.0, -0.24909197244573e-22],
    [1.0, 0.0, -0.26107636489332],
    [1.0, 1.0, 0.22592965981586],
    [1.0, 2.0, -0.64256463395226e-01],
    [1.0, 3.0, 0.78876289270526e-02],
    [1.0, 12.0, 0.35672110607366e-09],
    [1.0, 31.0, 0.17332496994895e-23],
    [2.0, 0.0, 0.56608900654837e-03],
    [2.0, 1.0, -0.32635483139717e-03],
    [2.0, 2.0, 0.44778286690632e-04],
    [2.0, 9.0, -0.51322156908507e-09],
    [2.0, 31.0, -0.42522657042207e-25],
    [3.0, 10.0, 0.26400441360689e-12],
    [3.0, 32.0, 0.78124600459723e-28],
    [4.0, 32.0, -0.30732199903668e-30],
];

pub mod gamma_dimensionless_specific_gibbs_free_energy;
pub use gamma_dimensionless_specific_gibbs_free_energy::*;

/// derivatives for dimensionless gibbs free energy
/// used to calculate specific volume, entropy,
/// internal energy, enthalpy,
/// cp, cv 
/// speed of sound 
/// isentropic exopnent
///
/// isobaric cubic exapnsion coeff 
/// isothermal compressibility
pub mod gamma_derivatives;
pub use gamma_derivatives::*;

/// intensive properties caluclated using the gamma_derivatives
///
/// these include 
/// specific volume 
/// specific enthalpy
/// specific internal energy 
/// specific entropy
/// specific cp 
/// specific cv 
/// speed of sound 
/// isentropic exponent (not done)
/// isobaric cubic expansion coeff (not done)
/// isothermal compressibility (not done)
pub mod intensive_properties;
pub use intensive_properties::*;


/// important tests to ensure things are working correctly
#[cfg(test)]
mod tests;

use uom::si::f64::*;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;


/// Returns the region-1 tau (dimensionless temperature)
/// Pressure is assumed to be in Pa
pub fn tau_1(t: ThermodynamicTemperature) -> f64 {
    // Temperature is assumed to be in K
    let t_kelvin = t.get::<kelvin>();
    1386.0 / t_kelvin
}

/// Returns the region-1 pi (dimensionless pressure)
/// Temperature is assumed to be in K
pub fn pi_1(p: Pressure) -> f64 {

    let p_pascals = p.get::<pascal>();
    // Pressure is assumed to be in Pa
    p_pascals / (16.53e6)
}



/// contains code and functions for backward equations for region 1 
/// pressure and enthalpy (p,h) flash
pub mod backward_eqn_ph_1;
pub use backward_eqn_ph_1::*;

/// contains code and functions for backwards equations 
/// pressure and entropy (p,s) flash 
pub mod backward_eqn_ps_1;
pub use backward_eqn_ps_1::*;

pub mod backward_eqn_hs_1;
pub use backward_eqn_hs_1::*;
