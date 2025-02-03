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

/// Returns the region-5 sound velocity
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn w_tp_5(t: ThermodynamicTemperature, p: Pressure) -> Velocity {
    let tau = tau_5(t);
    let pi = pi_5(p);
    let num = 1.0 + 2.0 * pi * gamma_pi_5_res(t, p) + pi.powi(2) * gamma_pi_5_res(t, p).powi(2);
    let subnum = (1.0 + pi * gamma_pi_5_res(t, p) - tau * pi * gamma_pi_tau_5_res(t, p)).powi(2);
    let subden = tau.powi(2) * (gamma_tau_tau_5_ideal(t, p) + gamma_tau_tau_5_res(t, p));
    let den = 1.0 - pi.powi(2) * gamma_pi_pi_5_res(t, p) + subnum / subden;
    ((specific_gas_constant_of_water()  * t) * num / den).sqrt()
}

/// Returns the region-5 isentropic exponent
pub fn kappa_tp_5(t: ThermodynamicTemperature, p: Pressure) -> Ratio {
    let tau = tau_5(t);
    let pi = pi_5(p);
    let num = 1.0 + 2.0 * pi * gamma_pi_5_res(t, p) + pi.powi(2) * gamma_pi_5_res(t, p).powi(2);
    let subnum = (1.0 + pi * gamma_pi_5_res(t, p) - tau * pi * gamma_pi_tau_5_res(t, p)).powi(2);
    let subden = tau.powi(2) * (gamma_tau_tau_5_ideal(t, p) + gamma_tau_tau_5_res(t, p));
    let den = (1.0 - pi.powi(2) * gamma_pi_pi_5_res(t, p) 
        + subnum / subden) * pi * (gamma_pi_5_ideal(t, p) + gamma_pi_5_res(t, p));

    return (num/den).into();
}


/// Returns the region-5 isobaric cubic expansion coeff
pub fn alpha_v_tp_5(t: ThermodynamicTemperature, p: Pressure) -> TemperatureCoefficient {
    let tau = tau_5(t);
    let pi = pi_5(p);
    let one_over_t: TemperatureCoefficient = 
        t.recip();
    let num = 1.0 + pi * gamma_pi_5_res(t, p) - tau * pi * gamma_pi_tau_5_res(t, p);
    let den = 1.0 + pi * gamma_pi_5_res(t, p);

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
/// Returns the region-5 isobaric isothermal compressibility
pub fn kappa_t_tp_5(t: ThermodynamicTemperature, p: Pressure) -> InversePressure {
    let pi = pi_5(p);
    let num = 1.0 - pi.powi(2) * gamma_pi_pi_5_res(t, p);
    let den = 1.0 + pi * gamma_pi_5_res(t, p);

    return (num/den)/p;

}
