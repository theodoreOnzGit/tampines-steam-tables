const REGION_5_COEFFS_RES: [[f64; 3]; 6] = [
    [1.0, 1.0, 0.15736404855259e-2],
    [1.0, 2.0, 0.90153761673944e-3],
    [1.0, 3.0, -0.50270077677648e-2],
    [2.0, 3.0, 0.22440037409485e-5],
    [2.0, 9.0, -0.41163275453471e-5],
    [3.0, 7.0, 0.37919454822955e-7],
];

const REGION_5_COEFFS_IDEAL: [[f64; 2]; 6] = [
    [0.0, -0.13179983674201e2],
    [1.0, 0.68540841634434e1],
    [-3.0, -0.24805148933466e-1],
    [-2.0, 0.36901534980333],
    [-1.0, -0.31161318213925e1],
    [2.0, -0.32961626538917],
];

pub mod dimensionless_tau_and_pi;
pub use dimensionless_tau_and_pi::*;

pub mod gamma_ideal_gas_plus_derivatives;
pub use gamma_ideal_gas_plus_derivatives::*;

pub mod gamma_residual_plus_derivatives;
pub use gamma_residual_plus_derivatives::*;

pub mod intensive_properties;
pub use intensive_properties::*;

#[cfg(test)]
mod tests;

