use uom::ConstZero;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::{TampinesSteamTableCV};

#[cfg(test)]
mod tests;


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

    // we should handle two edge cases, 
    // one outside the vicinity of choked flow,
    // one near it.
    //

    let dp_times_denominator = da_over_a * (rho1 * v1 * v1);
    let dv_times_denominator = -da_over_a * v1;

    // now we have singularity issue, near sonic flow, choked flow occurs
    // in such a regime, dv approaches zero.

    let m_value = mach_number_at_inlet.get::<ratio>();
    if m_value > 0.99 && m_value < 1.01 {
        // At throat: use isentropic relations, not dp/dA formula
        // The throat condition is: M = 1, which fixes the relationship
        //
        // dp and da go to zero here

        let speed_of_sound = state_1.get_speed_of_sound();
        let dv = speed_of_sound - v1;
        let critical_pressure_ratio: Ratio 
            = state_1.get_critical_pressure_ratio();

        let critical_pressure = critical_pressure_ratio * p1;


        let dp = critical_pressure - p1;

        return (dp,dv);

    }
    
    let dp = dp_times_denominator / denominator;
    let dv = dv_times_denominator / denominator;
    
    (dp, dv)
}

/// Given a step size, obtain the outlet pressure, velocity and 
/// enthalpy of the system using incremental quasi-1D analysis
///
/// da_by_a_max: Maximum fractional area change per step (e.g., 0.05 for 5%)
///              Should be small to maintain accuracy of differential relations
///
/// note: written originally, used AI to check
///
/// this is good for mach number 0-0.7 and 1.3 onwards
pub fn get_outlet_pressure_velocity_enthalpy_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity,
    da_by_a_max: Ratio) -> (Pressure, Velocity, AvailableEnergy) {

    // This process is isentropic, hence obtain the entropy
    let s1 = s_ph_eqm(p1, h1);

    // Total fractional area change needed
    let da_by_a_total: Ratio = (a2 - a1) / a1;  
    
    // Track remaining area change (use absolute value for comparison)
    let mut da_by_a_remaining = da_by_a_total.abs();
    
    // Determine sign of area change
    let area_change_sign: f64 = if da_by_a_total >= Ratio::ZERO { 
        1.0 
    } else { 
        -1.0 
    };

    // Set initial parameters 
    let mut area_before = a1;
    let mut p_intermediate = p1;
    let mut h_intermediate = h1;
    let mut v_intermediate = v1;

    // Tolerance for convergence
    let tolerance = Ratio::new::<ratio>(1e-6);

    while da_by_a_remaining > tolerance {
        
        // Determine step size for this iteration
        let da_by_a_step = if da_by_a_max < da_by_a_remaining {
            da_by_a_max
        } else {
            da_by_a_remaining
        };
        
        // Apply sign to step
        let da_by_a_signed = Ratio::new::<ratio>(
            da_by_a_step.get::<ratio>() * area_change_sign
        );
        
        // Calculate area after this step
        let area_after = area_before * (Ratio::new::<ratio>(1.0) + da_by_a_signed);

        // Get pressure and velocity changes
        let (dp, dv) = get_dp_dv_isentropic_nozzle_diffuser(
            area_before, 
            area_after, 
            p_intermediate, 
            h_intermediate, 
            v_intermediate
        );

        // Update state
        p_intermediate = p_intermediate + dp;
        v_intermediate = v_intermediate + dv;
        h_intermediate = h_ps_eqm(p_intermediate, s1);
        
        // Move to next step
        area_before = area_after;
        da_by_a_remaining = da_by_a_remaining - da_by_a_step;
    }

    let p2 = p_intermediate;
    let h2 = h_intermediate;
    let v2 = v_intermediate;

    (p2, v2, h2)
}


pub mod momentum_balance_rayleigh_line;
