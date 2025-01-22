use crate::constants::R_KJ_PER_KG_KELVIN;
use crate::region_2_saturated_water::{pi_2, tau_2};
use uom::si::{f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};
use super::{gamma_metastable_2_ideal, gamma_metastable_2_res, gamma_metastable_pi_2_ideal, gamma_metastable_pi_2_res, gamma_metastable_pi_pi_2_res, gamma_metastable_pi_tau_2_res, gamma_metastable_tau_2_ideal, gamma_metastable_tau_2_res, gamma_metastable_tau_tau_2_ideal, gamma_metastable_tau_tau_2_res};


#[inline]
pub fn specific_gas_constant_of_water() -> SpecificHeatCapacity {
    let r = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        R_KJ_PER_KG_KELVIN
    );

    r
}
/// Returns the region-2 specific volume
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn v_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    ((specific_gas_constant_of_water() ) * t / p) * pi_2(p) * (gamma_metastable_pi_2_ideal(t, p) + gamma_metastable_pi_2_res(t, p))
}

/// Returns the region-2 enthalpy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn h_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    specific_gas_constant_of_water() * t * tau_2(t) * (gamma_metastable_tau_2_ideal(t, p) + gamma_metastable_tau_2_res(t, p))
}

/// Returns the region-2 internal energy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn u_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    let tau = tau_2(t);
    let pi = pi_2(p);
    let tau_term = gamma_metastable_tau_2_ideal(t, p) + gamma_metastable_tau_2_res(t, p);
    let pi_term = gamma_metastable_pi_2_ideal(t, p) + gamma_metastable_pi_2_res(t, p);
    specific_gas_constant_of_water() * t * (tau * tau_term - pi * pi_term)
}

/// Returns the region-2 entropy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn s_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let tau = tau_2(t);
    let tau_term = gamma_metastable_tau_2_ideal(t, p) + gamma_metastable_tau_2_res(t, p);
    let pi_term = gamma_metastable_2_ideal(t, p) + gamma_metastable_2_res(t, p);
    specific_gas_constant_of_water() * (tau * tau_term - pi_term)
}

/// Returns the region-2 isobaric specific heat
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cp_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let tau = tau_2(t);
    -specific_gas_constant_of_water() * tau.powi(2) * (gamma_metastable_tau_tau_2_ideal(t, p) + gamma_metastable_tau_tau_2_res(t, p))
}

/// Returns the region-2 isochoric specific heat
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cv_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let tau = tau_2(t);
    let pi = pi_2(p);
    let num = (1.0 + pi * gamma_metastable_pi_2_res(t, p) - tau * pi * gamma_metastable_pi_tau_2_res(t, p)).powi(2);
    let den = 1.0 - pi.powi(2) * gamma_metastable_pi_pi_2_res(t, p);
    cp_tp_2_metastable(t, p) - specific_gas_constant_of_water() * num / den
}

/// Returns the region-2 sound velocity
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn w_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> Velocity {
    let tau = tau_2(t);
    let pi = pi_2(p);
    let num = 1.0 + 2.0 * pi * gamma_metastable_pi_2_res(t, p) + pi.powi(2) * gamma_metastable_pi_2_res(t, p).powi(2);
    let subnum = (1.0 + pi * gamma_metastable_pi_2_res(t, p) - tau * pi * gamma_metastable_pi_tau_2_res(t, p)).powi(2);
    let subden = tau.powi(2) * (gamma_metastable_tau_tau_2_ideal(t, p) + gamma_metastable_tau_tau_2_res(t, p));
    let den = 1.0 - pi.powi(2) * gamma_metastable_pi_pi_2_res(t, p) + subnum / subden;
    ((specific_gas_constant_of_water()  * t) * num / den).sqrt()
}

/// Returns the region-2 isentropic exponent
pub fn kappa_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> Ratio {
    let tau = tau_2(t);
    let pi = pi_2(p);
    let num = 1.0 + 2.0 * pi * gamma_metastable_pi_2_res(t, p) + pi.powi(2) * gamma_metastable_pi_2_res(t, p).powi(2);
    let subnum = (1.0 + pi * gamma_metastable_pi_2_res(t, p) - tau * pi * gamma_metastable_pi_tau_2_res(t, p)).powi(2);
    let subden = tau.powi(2) * (gamma_metastable_tau_tau_2_ideal(t, p) + gamma_metastable_tau_tau_2_res(t, p));
    let den = (1.0 - pi.powi(2) * gamma_metastable_pi_pi_2_res(t, p) 
        + subnum / subden) * pi * (gamma_metastable_pi_2_ideal(t, p) + gamma_metastable_pi_2_res(t, p));

    return (num/den).into();
}


/// Returns the region-2 isobaric cubic expansion coeff
pub fn alpha_v_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> TemperatureCoefficient {
    let tau = tau_2(t);
    let pi = pi_2(p);
    let one_over_t: TemperatureCoefficient = 
        t.recip();
    let num = 1.0 + pi * gamma_metastable_pi_2_res(t, p) - tau * pi * gamma_metastable_pi_tau_2_res(t, p);
    let den = 1.0 + pi * gamma_metastable_pi_2_res(t, p);

    return one_over_t * num/den;

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
/// Returns the region-1 isobaric isothermal compressibility
pub fn kappa_t_tp_2_metastable(t: ThermodynamicTemperature, p: Pressure) -> InversePressure {
    let pi = pi_2(p);
    let num = 1.0 - pi.powi(2) * gamma_metastable_pi_pi_2_res(t, p);
    let den = 1.0 + pi * gamma_metastable_pi_2_res(t, p);

    return (num/den)/p;

}

