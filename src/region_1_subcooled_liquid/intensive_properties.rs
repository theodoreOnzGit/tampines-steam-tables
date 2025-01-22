use uom::si::temperature_coefficient::per_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::{f64::*, ratio::ratio};

use crate::constants::specific_gas_constant_of_water;

use super::{gamma_1, gamma_pi_1, gamma_pi_pi_1, gamma_pi_tau_1, gamma_tau_1, gamma_tau_tau_1, pi_1, tau_1};


/// Returns the region-1 specific enthalpy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn h_tp_1(t: ThermodynamicTemperature, p: Pressure) -> 
AvailableEnergy {
    specific_gas_constant_of_water() * t * tau_1(t) * gamma_tau_1(t, p)
}

/// Returns the region-1 specific volume
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn v_tp_1(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    // in rust_steam 
    // The multiplication by 1000 is necessary to convert R from kJ/kg.K to J/kg.K
    // but the uom package takes care of that so we are not dealing with this anymore
    ((specific_gas_constant_of_water() ) * t / p) * pi_1(p) * gamma_pi_1(t, p)
}

/// Returns the region-1 specific internal energy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn u_tp_1(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    specific_gas_constant_of_water() * t * (tau_1(t) * gamma_tau_1(t, p) - pi_1(p) * gamma_pi_1(t, p))
}

/// Returns the region-1 specific entropy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
///
/// units are same as cp
pub fn s_tp_1(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    specific_gas_constant_of_water() * (tau_1(t) * gamma_tau_1(t, p) - gamma_1(t, p))
}

/// Returns the region-1 specific isobaric heat capacity
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cp_tp_1(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let tau = tau_1(t);
    specific_gas_constant_of_water() * (-tau.powi(2) * gamma_tau_tau_1(t, p))
}

/// Returns the region-1 specific isochoric heat capacity
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cv_tp_1(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let tau = tau_1(t);
    let corr = (gamma_pi_1(t, p) - tau * gamma_pi_tau_1(t, p)).powi(2) / gamma_pi_pi_1(t, p);
    specific_gas_constant_of_water() * (-tau.powi(2) * gamma_tau_tau_1(t, p) + corr)
}

/// Returns the region-1 speed of sound
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn w_tp_1(t: ThermodynamicTemperature, p: Pressure) -> Velocity {
    let tau = tau_1(t);
    let gamma_pi = gamma_pi_1(t, p);
    let gamma_pi_tau = gamma_pi_tau_1(t, p);
    let gamma_pi_pi = gamma_pi_pi_1(t, p);
    let gamma_tau_tau = gamma_tau_tau_1(t, p);
    let term = (gamma_pi - tau * gamma_pi_tau).powi(2) / (tau.powi(2) * gamma_tau_tau);

    // in rust_steam 
    // The multiplication by 1000 is necessary to convert R from kJ/kg.K to J/kg.K
    // however, the units of measure crate takes care of it
    let square = (specific_gas_constant_of_water() ) * t * (gamma_pi.powi(2) / (term - gamma_pi_pi));
    square.sqrt()
}


/// Returns the region-1 isentropic exponent
pub fn kappa_tp_1(t: ThermodynamicTemperature, p: Pressure) -> Ratio {
    let tau = tau_1(t);
    let pi = pi_1(p);
    let gamma_pi = gamma_pi_1(t, p);
    let gamma_pi_tau = gamma_pi_tau_1(t, p);
    let gamma_pi_pi = gamma_pi_pi_1(t, p);
    let gamma_tau_tau = gamma_tau_tau_1(t, p);
    let denominator = (gamma_pi - tau * gamma_pi_tau).powi(2) / (tau.powi(2) * gamma_tau_tau)*pi - pi * gamma_pi_pi;

    let numerator = gamma_pi;

    return (numerator/denominator).into();

}


/// Returns the region-1 isobaric cubic expansion coeff
pub fn alpha_v_tp_1(t: ThermodynamicTemperature, p: Pressure) -> TemperatureCoefficient {
    let tau = tau_1(t);
    let gamma_pi = gamma_pi_1(t, p);
    let gamma_pi_tau = gamma_pi_tau_1(t, p);

    let dimensionless_alpha: Ratio = Ratio::new::<ratio>(1.0 - tau * gamma_pi_tau / gamma_pi);
    let t_kelvin = t.get::<kelvin>();

    return dimensionless_alpha * TemperatureCoefficient::new::<per_kelvin>(t_kelvin.recip());

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
pub fn kappa_t_tp_1(t: ThermodynamicTemperature, p: Pressure) -> InversePressure {
    let pi = pi_1(p);
    let gamma_pi = gamma_pi_1(t, p);
    let gamma_pi_pi = gamma_pi_pi_1(t, p);

    let dimensionless_kappa_t: Ratio = -Ratio::new::<ratio>( pi * gamma_pi_pi / gamma_pi );

    return dimensionless_kappa_t/p;

}
