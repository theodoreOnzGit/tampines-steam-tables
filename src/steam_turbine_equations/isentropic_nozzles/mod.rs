use uom::ConstZero;
use uom::si::f64::*;
use uom::si::force::newton;
use uom::si::pressure::{bar, millibar, pascal};
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::hs_flash_eqm::v_hs_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::{s_ph_eqm, w_ph_eqm};
use crate::prelude::functional_programming::ps_flash_eqm::{h_ps_eqm, v_ps_eqm};
use crate::prelude::{TampinesSteamTableCV};


// for isentropic/adiabatic diffusers/nozzles, we follow the mach number
// formulae 
//
pub fn get_dp_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity) -> Pressure {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let rho1: MassDensity = state_1.get_rho();
    let mach_number_at_inlet: Ratio = state_1.get_mach_number(v1);
    let da = a2 - a1;
    let ratio_one = Ratio::new::<ratio>(1.0);

    let dp = da/a1 * (rho1 * v1 * v1)/(ratio_one - mach_number_at_inlet *mach_number_at_inlet);

    dp

}


// for isentropic/adiabatic diffusers/nozzles, we follow the mach number
// formulae 
//
pub fn get_dv_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity) -> Velocity {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let mach_number_at_inlet: Ratio = state_1.get_mach_number(v1);
    let da = a2 - a1;
    let ratio_one = Ratio::new::<ratio>(1.0);

    let dv = -da/a1 * v1/(ratio_one - mach_number_at_inlet *mach_number_at_inlet);

    dv

}

/// Returns both pressure and velocity changes for isentropic nozzle/diffuser
/// More efficient than calling get_dp and get_dv separately
pub fn get_dp_dv_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity) -> (Pressure, Velocity) {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let rho1: MassDensity = state_1.get_rho();
    let mach_number_at_inlet: Ratio = state_1.get_mach_number(v1);
    let m_squared = mach_number_at_inlet * mach_number_at_inlet;
    let ratio_one = Ratio::new::<ratio>(1.0);
    
    let da_over_a = (a2 - a1) / a1;
    let denominator = ratio_one - m_squared;
    
    let dp = da_over_a * (rho1 * v1 * v1) / denominator;
    let dv = -da_over_a * v1 / denominator;
    
    (dp, dv)
}

/// given a step size, obtain the outlet pressure, velocity and 
/// enthalpy of the system
pub fn get_outlet_pressure_velocity_enthalpy_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity,
    da_by_a_max: Ratio) -> (Pressure, Velocity, AvailableEnergy) {

    // this process is isentropic, hence obtain the entropy
    let s1 = s_ph_eqm(p1, h1);

    let da_by_a_total: Ratio = (a1-a2)/a1;

    let mut da_by_a_remaining = da_by_a_total;

    // set some parameters 
    let mut area_before = a1;
    let mut da_by_a = da_by_a_max;


    if da_by_a >= da_by_a_remaining {
        da_by_a = da_by_a_remaining;
    }

    let mut area_after = area_before + area_before * da_by_a_max;
    let mut p_intermediate = p1;
    let mut h_intermediate = h1;
    let mut v_intermediate = v1;

    while da_by_a_remaining >= Ratio::ZERO {

        let (dp, dv) = get_dp_dv_isentropic_nozzle_diffuser(
            area_before, area_after, p_intermediate, h_intermediate, v_intermediate
        );

        // reset the parameters 
        p_intermediate += dp;
        v_intermediate += dv;

        if da_by_a >= da_by_a_remaining {
            da_by_a = da_by_a_remaining;
        }

        area_before = area_after;
        area_after = area_after + area_after * da_by_a;
        h_intermediate = h_ps_eqm(p_intermediate, s1);

        da_by_a_remaining -= da_by_a;
    }

    let p2 = p_intermediate;
    let h2 = h_intermediate;
    let v2 = v_intermediate;

    return (p2, v2, h2);
}


