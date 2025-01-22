use super::{region_4_coeff_index, theta_dimensionless_temp_4};
use uom::si::{f64::*, pressure::megapascal};

/// returns sat pressure in region 4
pub fn sat_pressure_4(t: ThermodynamicTemperature) -> Pressure {
    let theta = theta_dimensionless_temp_4(t);
    let ref_p = Pressure::new::<megapascal>(1.0);

    return ref_p * dimensionless_sat_pressure(theta);

}

#[inline]
fn dimensionless_sat_pressure(theta: f64) -> f64 {
    let a = coeff_a(theta);
    let b = coeff_b(theta);
    let c = coeff_c(theta);

    let num = 2.0 * c;
    let den = -b + (b.powi(2) - 4.0 * a * c).sqrt();

    return (num/den).powi(4);
}

// in sat pressure eqn, it's A
#[inline]
fn coeff_a(theta: f64) -> f64 {
    let n1 = region_4_coeff_index(1);
    let n2 = region_4_coeff_index(2);

    theta.powi(2) + n1 * theta + n2

}

// in sat pressure eqn, it's B
#[inline]
fn coeff_b(theta: f64) -> f64 {
    let n3 = region_4_coeff_index(3);
    let n4 = region_4_coeff_index(4);
    let n5 = region_4_coeff_index(5);

    n3 * theta.powi(2) + n4 * theta + n5

}

// in sat pressure eqn, it's C
#[inline]
fn coeff_c(theta: f64) -> f64 {
    let n6 = region_4_coeff_index(6);
    let n7 = region_4_coeff_index(7);
    let n8 = region_4_coeff_index(8);

    n6 * theta.powi(2) + n7 * theta + n8

}
