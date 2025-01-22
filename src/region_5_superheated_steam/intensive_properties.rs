use crate::constants::specific_gas_constant_of_water;
use uom::si::f64::*;

use super::{gamma_5_ideal, gamma_5_res, gamma_pi_5_ideal, gamma_pi_5_res, gamma_pi_pi_5_res, gamma_pi_tau_5_res, gamma_tau_5_ideal, gamma_tau_5_res, gamma_tau_tau_5_ideal, gamma_tau_tau_5_res, pi_5, tau_5};
/// Returns the region-5 specific volume
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn v_tp_5(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    ((specific_gas_constant_of_water()) * t / p) * pi_5(p) * (gamma_pi_5_ideal(t, p) + gamma_pi_5_res(t, p))
}

/// Returns the region-5 enthalpy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn h_tp_5(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    specific_gas_constant_of_water() * t * tau_5(t) * (gamma_tau_5_ideal(t, p) + gamma_tau_5_res(t, p))
}

/// Returns the region-5 internal energy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn u_tp_5(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    let tau: f64 = tau_5(t);
    let pi: f64 = pi_5(p);
    specific_gas_constant_of_water()
        * t
        * (tau * (gamma_tau_5_ideal(t, p) + gamma_tau_5_res(t, p))
            - pi * (gamma_pi_5_ideal(t, p) + gamma_pi_5_res(t, p)))
}

/// Returns the region-5 entropy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn s_tp_5(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let tau = tau_5(t);
    specific_gas_constant_of_water()
        * (tau * (gamma_tau_5_ideal(t, p) + gamma_tau_5_res(t, p))
            - (gamma_5_ideal(t, p) + gamma_5_res(t, p)))
}

/// Returns the region-5 isobaric specific heat
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cp_tp_5(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    -specific_gas_constant_of_water() * tau_5(t).powi(2) * (gamma_tau_tau_5_ideal(t, p) + gamma_tau_tau_5_res(t, p))
}

/// Returns the region-5 isochoric specific heat
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cv_tp_5(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let pi: f64 = pi_5(p);
    cp_tp_5(t, p)
        - specific_gas_constant_of_water()
            * (((1.0 + pi * gamma_pi_5_res(t, p) - tau_5(t) * pi * gamma_pi_tau_5_res(t, p))
                .powi(2))
                / (1.0 - pi.powi(2) * gamma_pi_pi_5_res(t, p)))
}

