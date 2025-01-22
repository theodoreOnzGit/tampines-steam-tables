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

/// speed of sound in region 3
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


// to make the inverse pressure type 
// it is m s^2 / kg 
use uom::si::{ISQ, SI, Quantity};
use uom::typenum::{Z0, P1, P2, N1};

// quantity is defined
// ## Generic Parameters
// * `L`: Length dimension.
// * `M`: Mass dimension.
// * `T`: Time dimension.
// * `I`: Electric current dimension.
// * `Th`: Thermodynamic temperature dimension.
// * `N`: Amount of substance dimension.
// * `J`: Luminous intensity dimension.
// * `K`: Kind.
pub type InversePressure = Quantity<ISQ<P1, N1, P2, Z0, Z0, Z0, Z0>, SI<f64>, f64>;

/// isentropic exponent in region 3 
pub fn kappa_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {

    let delta = delta_3(rho);
    let tau = tau_3(t);
    let phi_delta = phi_delta_3(rho, t);
    let phi_delta_delta = phi_delta_delta_3(rho, t);
    let phi_delta_tau = phi_delta_tau_3(rho, t);
    let phi_tau_tau = phi_tau_tau_3(rho, t);

    let first_term = 2.0 + delta * phi_delta_delta/phi_delta;

    let second_term_num = (delta * phi_delta - delta * tau * phi_delta_tau).powi(2);
    let second_term_den = delta * tau.powi(2) * phi_delta * phi_tau_tau;

    return first_term - second_term_num/second_term_den;

}


/// Returns the region-3 isobaric cubic expansion coeff
pub fn alpha_v_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> TemperatureCoefficient {

    let delta = delta_3(rho);
    let tau = tau_3(t);
    let phi_delta = phi_delta_3(rho, t);
    let phi_delta_tau = phi_delta_tau_3(rho, t);
    let phi_tau_tau = phi_tau_tau_3(rho, t);

    let num = phi_delta - tau * phi_delta_tau;
    let den = 2.0 * phi_delta + delta * phi_tau_tau;

    return t.recip() * (num/den);

}

/// Returns the region-3 isothermal compressibility
pub fn kappa_t_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> InversePressure {

    let r = specific_gas_constant_of_water();
    let rho_r_t: Pressure = rho * r * t;

    let delta = delta_3(rho);
    let phi_delta = phi_delta_3(rho, t);
    let phi_delta_delta = phi_delta_delta_3(rho, t);

    let den = 2.0 * delta * phi_delta + delta.powi(2) * phi_delta_delta;

    return rho_r_t.recip() * den.recip();


}

/// Returns the region-3 relative pressure coefficient
pub fn alpha_p_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> TemperatureCoefficient {

    let tau = tau_3(t);
    let phi_delta = phi_delta_3(rho, t);
    let phi_delta_tau = phi_delta_tau_3(rho, t);

    return t.recip() * (1.0 - tau * phi_delta_tau/phi_delta );

}


/// Returns the region-3 isothermal stress coefficient
pub fn beta_p_rho_t_3(rho: MassDensity, t: ThermodynamicTemperature) -> MassDensity {

    let delta = delta_3(rho);
    let phi_delta = phi_delta_3(rho, t);
    let phi_delta_delta = phi_delta_delta_3(rho, t);

    let num = 2.0 + delta * phi_delta_delta / phi_delta;

    return rho * num;

}
