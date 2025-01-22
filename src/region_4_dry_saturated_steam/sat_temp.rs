use super::{beta_dimensionless_pressure_4, region_4_coeff_index};
use uom::si::{f64::*, thermodynamic_temperature::kelvin};

pub fn sat_temp_4(p: Pressure) -> ThermodynamicTemperature {

    let beta = beta_dimensionless_pressure_4(p);
    let ref_t = ThermodynamicTemperature::new::<kelvin>(1.0);
    return dimensionless_sat_temp(beta) * ref_t;

}

// dimensionless sat temp
#[inline]
fn dimensionless_sat_temp(beta: f64) -> f64 {
    let d = coeff_d(beta);

    let n9 = region_4_coeff_index(9);
    let n10 = region_4_coeff_index(10);

    let num = n10 + d - ( (n10 + d).powi(2) - 4.0 * (n9 + n10 * d) ).sqrt();

    let _den = 2.0;

    // normally just numerator over denominator, but just multiply by 0.5
    return num * 0.5;


}


// in sat temp eqn, it's D
#[inline]
fn coeff_d(beta: f64) -> f64 {
    let e = coeff_e(beta);
    let f = coeff_f(beta);
    let g = coeff_g(beta);

    let num = 2.0 * g;
    let den = -f - (f.powi(2) - 4.0 * e * g).sqrt();
    return num/den;

}

// in sat temp eqn, it's E
#[inline]
fn coeff_e(beta: f64) -> f64 {
    let n3 = region_4_coeff_index(3);
    let n6 = region_4_coeff_index(6);

    beta.powi(2) + n3 * beta + n6

}

// in sat temp eqn, it's F
#[inline]
fn coeff_f(beta: f64) -> f64 {
    let n1 = region_4_coeff_index(1);
    let n4 = region_4_coeff_index(4);
    let n7 = region_4_coeff_index(7);

    n1 * beta.powi(2) + n4 * beta + n7

}

// in sat temp eqn, it's G
#[inline]
fn coeff_g(beta: f64) -> f64 {
    let n2 = region_4_coeff_index(2);
    let n5 = region_4_coeff_index(5);
    let n8 = region_4_coeff_index(8);

    n2 * beta.powi(2) + n5 * beta + n8

}

