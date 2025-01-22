use crate::constants::R_KJ_PER_KG_KELVIN;

use super::{delta_3, phi_3, phi_delta_3, phi_delta_delta_3, phi_delta_tau_3, phi_tau_3, phi_tau_tau_3, tau_3};
use uom::si::{f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};
#[inline]
pub fn specific_gas_constant_of_water() -> SpecificHeatCapacity {
    let r = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        R_KJ_PER_KG_KELVIN
    );

    r
}


/// Returns the pressure given t and rho
/// Temperature is assumed to be in K
/// density is assumed to be in kg/m^3
pub fn p_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> Pressure {
    rho * (specific_gas_constant_of_water()) * t * delta_3(rho) * phi_delta_3(rho, t)
}

/// Returns the internal energy given t and rho
/// Temperature is assumed to be in K
/// density is assumed to be in kg/m^3
pub fn u_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> AvailableEnergy {
    tau_3(t) * phi_tau_3(rho, t) * specific_gas_constant_of_water() * t
}

/// Returns the entropy given t and rho
/// Temperature is assumed to be in K
/// density is assumed to be in kg/m^3
pub fn s_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> SpecificHeatCapacity {
    (tau_3(t) * phi_tau_3(rho, t) - phi_3(rho, t)) * specific_gas_constant_of_water()
}

/// Returns the enthalpy given t and rho
/// Temperature is assumed to be in K
/// density is assumed to be in kg/m^3
pub fn h_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> AvailableEnergy {
    (tau_3(t) * phi_tau_3(rho, t) + delta_3(rho) * phi_delta_3(rho, t)) * specific_gas_constant_of_water() * t
}

/// Returns the isochoric specific heat given t and rho
/// Temperature is assumed to be in K
/// density is assumed to be in kg/m^3
pub fn cv_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> SpecificHeatCapacity {
    -tau_3(t).powi(2) * phi_tau_tau_3(rho, t) * specific_gas_constant_of_water()
}

/// Returns the isobaric specific heat given t and rho
/// Temperature is assumed to be in K
/// density is assumed to be in kg/m^3
pub fn cp_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> SpecificHeatCapacity {
    (-tau_3(t).powi(2) * phi_tau_tau_3(rho, t)
        + ((delta_3(rho) * phi_delta_3(rho, t)
            - delta_3(rho) * tau_3(t) * phi_delta_tau_3(rho, t))
        .powi(2)
            / (2.0 * delta_3(rho) * phi_delta_3(rho, t)
                + delta_3(rho).powi(2) * phi_delta_delta_3(rho, t))))
        * specific_gas_constant_of_water()
}

pub fn w_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> Velocity {
    ((2.0 * delta_3(rho) * phi_delta_3(rho, t) + delta_3(rho).powi(2) * phi_delta_delta_3(rho, t)
        - ((delta_3(rho) * phi_delta_3(rho, t)
            - delta_3(rho) * tau_3(t) * phi_delta_tau_3(rho, t))
        .powi(2)
            / (tau_3(t).powi(2) * phi_tau_tau_3(rho, t))))
        * specific_gas_constant_of_water()
        * t)
        .sqrt()
}
