use uom::si::{dynamic_viscosity::pascal_second, f64::*, ratio::ratio};

use crate::{constants::{rho_crit_water, t_crit_water}, interfaces::functional_programming::{ph_flash_eqm::{t_ph_eqm, v_ph_eqm}, pt_flash_eqm::{v_tp_eqm_single_phase, v_tp_eqm_two_phase}}};


const PSI_0_COEFFS: [[f64; 2]; 4] = [
    [1.0,  0.167_752e-1],
    [2.0,  0.220_462e-1],
    [3.0,  0.636_656_4e-2],
    [4.0,  -0.241_605e-2],
];

const PSI_1_COEFFS: [[f64; 4]; 21] = [
    [1.0, 0.0, 0.0, 0.520_094],
    [2.0, 0.0, 1.0, 0.850_895e-1],
    [3.0, 0.0, 2.0, -0.108_374e1],
    [4.0, 0.0, 3.0, -0.289_555],
    [5.0, 1.0, 0.0, 0.222_531],
    [6.0, 1.0, 1.0, 0.999_115],
    [7.0, 1.0, 2.0, 0.188_797e1],
    [8.0, 1.0, 3.0, 0.126_613e1],
    [9.0, 1.0, 5.0, 0.120_573],
    [10.0, 2.0, 0.0, -0.281_378],
    [11.0, 2.0, 1.0, -0.906_851],
    [12.0, 2.0, 2.0, -0.772_479],
    [13.0, 2.0, 3.0, -0.489_837],
    [14.0, 2.0, 4.0, -0.257_040],
    [15.0, 3.0, 0.0, 0.161_913],
    [16.0, 3.0, 1.0, 0.257_399],
    [17.0, 4.0, 0.0, -0.325_372e-1],
    [18.0, 4.0, 3.0, 0.698_452e-1],
    [19.0, 5.0, 4.0, 0.872_102e-2],
    [20.0, 6.0, 3.0, -0.435_673e-2],
    [21.0, 6.0, 5.0, -0.593_264e-3],
];

/// for viscosity estimates in two phase region
/// and single phase region
pub fn mu_tp_eqm_two_phase(t: ThermodynamicTemperature,
    p: Pressure,
    x: f64) -> DynamicViscosity {
    let rho = v_tp_eqm_two_phase(t, p, x).recip();
    let psi = psi_0_viscosity(t) * psi_1_viscosity(t, rho);
    let eta_star = DynamicViscosity::new::<pascal_second>(1.0e-6);

    return psi * eta_star;

}
/// for viscosity estimates in single phase region
pub fn mu_tp_eqm_single_phase(t: ThermodynamicTemperature,
    p: Pressure) -> DynamicViscosity {
    let rho = v_tp_eqm_single_phase(t, p).recip();
    let psi = psi_0_viscosity(t) * psi_1_viscosity(t, rho);
    let eta_star = DynamicViscosity::new::<pascal_second>(1.0e-6);

    return psi * eta_star;

}
/// for viscosity estimates in two phase region
/// and single phase region
pub fn mu_rho_t_eqm(t: ThermodynamicTemperature,
    rho: MassDensity,) -> DynamicViscosity {

    let psi = psi_0_viscosity(t) * psi_1_viscosity(t, rho);
    let eta_star = DynamicViscosity::new::<pascal_second>(1.0e-6);

    return psi * eta_star;

}
/// for viscosity estimates in two phase region
/// and single phase region
/// using enthalpy and pressure
pub fn mu_ph_eqm(p: Pressure,
    h: AvailableEnergy) -> DynamicViscosity {

    let t = t_ph_eqm(p, h);
    let rho = v_ph_eqm(p, h).recip();

    return mu_rho_t_eqm(t, rho);

}

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
pub(crate) fn psi_1_viscosity(
    t: ThermodynamicTemperature,
    rho: MassDensity) -> f64 {
    let t_crit_water = t_crit_water();
    let theta = t/t_crit_water;
    let theta_f64 = theta.get::<ratio>();

    let rho_crit_water = rho_crit_water();
    let delta = rho/rho_crit_water;
    let delta_f64 = delta.get::<ratio>();

    let mut exponent: f64 = 0.0;

    for coeffs in PSI_1_COEFFS {
        let _i = coeffs[0];
        let ii = coeffs[1];
        let ji = coeffs[2];
        let ni = coeffs[3];

        exponent += ni * (delta_f64 - 1.0).powf(ii) * (theta_f64.recip() - 1.0).powf(ji);

    }

    exponent *= delta_f64;

    return exponent.exp();


}

#[cfg(test)]
mod tests;
