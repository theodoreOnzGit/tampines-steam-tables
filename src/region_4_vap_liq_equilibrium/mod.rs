// note that i add underscores here to make it easier 
// to read
const REGION_4_COEFFS: [f64; 10] = [
    0.116_705_214_527_67e4,
    -0.724_213_167_032_06e6,
    -0.170_738_469_400_92e2,
    0.120_208_247_024_70e5,
    -0.323_255_503_223_33e7,
    0.149_151_086_135_30e2,
    -0.482_326_573_615_91e4,
    0.405_113_405_420_57e6,
    -0.238_555_575_678_49,
    0.650_175_348_447_98e3,
];
/// helps to index region 4 coefficients
/// because of the indexing things (start at 0 vs start at 1)
#[inline]
pub(crate) fn region_4_coeff_index(i: usize) -> f64{
    REGION_4_COEFFS[i-1]
}

/// saturation temperature equation 
pub mod sat_temp;
pub use sat_temp::*;

/// saturation pressure equation 
pub mod sat_pressure;
pub use sat_pressure::*;

use uom::si::{f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

/// returns dimensionless pressure 
/// (there is an exponent to the power of 1/4)
/// in region 4
pub fn beta_dimensionless_pressure_4(p: Pressure) -> f64 {

    let pressure_mpa = 1.0;
    let ref_p = Pressure::new::<megapascal>(pressure_mpa);

    let pressure_ratio: f64 = (p/ref_p).into();

    return pressure_ratio.powf(0.25);
    

}

/// returns dimensionless temp for region 4
pub fn theta_dimensionless_temp_4(t: ThermodynamicTemperature) -> f64 {
    let ref_t = ThermodynamicTemperature::new::<kelvin>(1.0);
    let temp_ratio: f64 = (t/ref_t).into();

    let n9 = region_4_coeff_index(9);
    let n10 = region_4_coeff_index(10);

    return temp_ratio + n9/(temp_ratio - n10);
}

/// backward equation T_s(h,s)
pub mod backward_eqn_hs_4;
pub use backward_eqn_hs_4::*;

// tests 
#[cfg(test)]
mod tests;
