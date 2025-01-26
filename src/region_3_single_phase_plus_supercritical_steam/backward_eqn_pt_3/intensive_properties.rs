use uom::si::{f64::*, ratio::ratio};

use crate::region_3_single_phase_plus_supercritical_steam::{alpha_v_rho_t_3, cp_rho_t_3, cv_rho_t_3, h_rho_t_3, kappa_rho_t_3, kappa_t_rho_t_3, s_rho_t_3, u_rho_t_3, w_rho_t_3, InversePressure};

use super::v_tp_3;

/// Returns the region-3 enthalpy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn h_tp_3(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    let v = v_tp_3(t, p);
    let rho = v.recip();

    h_rho_t_3(rho, t)
}

/// Returns the region-3 internal energy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn u_tp_3(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {

    let v = v_tp_3(t, p);
    let rho = v.recip();

    u_rho_t_3(rho, t)
}

/// Returns the region-3 entropy
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn s_tp_3(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let v = v_tp_3(t, p);
    let rho = v.recip();

    s_rho_t_3(rho, t)

}

/// Returns the region-3 isobaric specific heat
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cp_tp_3(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let v = v_tp_3(t, p);
    let rho = v.recip();

    cp_rho_t_3(rho, t)
}

/// Returns the region-3 isochoric specific heat
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn cv_tp_3(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {

    let v = v_tp_3(t, p);
    let rho = v.recip();
    cv_rho_t_3(rho, t)
}

/// Returns the region-3 sound velocity
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn w_tp_3(t: ThermodynamicTemperature, p: Pressure) -> Velocity {
    let v = v_tp_3(t, p);
    let rho = v.recip();
    w_rho_t_3(rho, t)

}

/// Returns the region-3 isentropic exponent
pub fn kappa_tp_3(t: ThermodynamicTemperature, p: Pressure) -> Ratio {

    let v = v_tp_3(t, p);
    let rho = v.recip();
    Ratio::new::<ratio>(kappa_rho_t_3(rho, t))
}


/// Returns the region-3 isobaric cubic expansion coeff
pub fn alpha_v_tp_3(t: ThermodynamicTemperature, p: Pressure) -> TemperatureCoefficient {
    let v = v_tp_3(t, p);
    let rho = v.recip();
    alpha_v_rho_t_3(rho, t)

}


/// Returns the region-3 isobaric isothermal compressibility
pub fn kappa_t_tp_3(t: ThermodynamicTemperature, p: Pressure) -> InversePressure {

    let v = v_tp_3(t, p);
    let rho = v.recip();
    kappa_t_rho_t_3(rho, t)
}
