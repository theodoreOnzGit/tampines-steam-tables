use uom::ConstZero;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::{TampinesSteamTableCV};

#[cfg(test)]
mod tests;

// so here is the problem statement. 
//
//
// I have a converging diverging nozzle 
//
// there is a SET 
// inlet area, 
// throat area,
// outlet area,
//
// the conditions at inlet are known (this could be mass flowrate and 
// inlet velocity and such)
//
// The pressure at the turbine outlet will be known 
// (condenser pressure)
//
// It is very likely either 
//
// 1. saturated steam
// 2. superheated steam
//
// Water can in theory flow through the tubrine but that is highly 
// unlikely in steam turbine
//
// The turbine will have inlet pressure, and outlet pressure at condenser
//
// There will be a series of nozzles and turbine blades
// before then arriving at the outlet
// 
// The steam turbine can be a stator rotor pair interlaced with 
// many small control volumes in between each stator and rotor
//
// 
// for such a case, then there will be an inlet control volume with 
// known state, 
//
// then a stator (nozzle) and rotor pair
//
// For such a case, this is a pressure driven flow 
//
// the pressure difference between inlet and outlet will drive the flow
// It will be treated like a flow between two pressure points.
// Mass flowrate will be such that the pressures will equalise.
//
// For now, an impulse turbine stage can look like this 
// take two adjacent control volumes ,
//
//
// obtain pressure difference between them 
// knowing outlet pressure, I would then take this to be the nozzle 
// outlet pressure, 
//
// 1. I would then solve for the nozzle flow using the pressure differences
// 2. for the turbine blade, there is torque exerted, and the velocity 
// will be slowed down on the way out. This is based on the velocity 
// coefficient of the impulse turbine.
// 3. the outlet stagnation enthalpy is then calculated to be the state 
// of the steam entering the next control volume
//
// when pressure differences are given, we solve for mass flowrate 
// and velocity individually.
//
// I suppose we could do finite differencing. 
// In doing so, we may run into numerical errors such as Courant number.
// May be needed to do a local timestepping of sorts.
// Do we use rhoPimpleFoam algorithm? I myself am not sure.
//
// But the equations would be sort of an exponential type yeah... 
//
// not sure if I want to go down this rabbit hole 
// but nevermind.
//
// Anyway for a nozzle, with pressure differences, and no specified 
// mass flowrate 
//
// 1. assume stagnation in the inlet. (neglect KE in the inlet), otherwise 
// add it in
// 2. calculate throat velocity (possibly sonic), does this work?
// 3. generally speaking, nozzles just obey mass and energy balance.
// Energy balance is:
// v_throat = (2.0 * (h0 - h_throat).sqrt()
//
// assuming inlet stagnation states are known (p0, h0), we can 
// calculate entropy. This will also give us the inlet density of 
// steam.
//
// Nozzle is isentropic. (s1 = s0 = s_throat)
//
// area of entrance and throat are known 
//
// rho1 * a1 * v1 = rho2 * a2 * v2
//
// from this, we can obtain:
//
// rho2 * v2/v1 = rho1 * a1 / (a2);
//
// rho2 * v2/v1 is constant from mass balance
//
// What we can do, is to guess a mass flowrate.
//
// After guessing mass flowrate, we obtain v1 first 
// now we have the product of rho2 * v2.
//
// We have to guess h_throat iteratively. we shall have nested loop (oops)
//
// Otherwise, guess h_throat iteratively, the limits will be between 
// v_throat is supersonic and v_throat is 0
//
// from then we can guess v_throat. 
//
// Better yet, just guess v_throat (v2) first.
//
// This makes it easy to guess h_throat. 
// When h_throat is known, s_throat is known, thermodynamic state is fully 
// known. 
//
// we then get a formula for v2/v1. 
//
// From this, we get v1. Which is the mass flowrate.
//
// Now, when v_throat is guessed, the pressure at the throat is also 
// guessed.
//
// We have to see whether the pressure at the throat is good for the 
// pressure bounds given.
//
// Once having a throat pressure, we can obtain an equation to get 
// the mass flowrate in the outside part of the nozzle given the 
// nozzle diverging part.
//
// When the mass flowrate of the nozzle and supersonic diffuser match,
// given the pressure difference set, 
// then we have a mass flowrate solution in the nozzle.
// 
//
//
//
//
//
// 
//
//
//
//
//
//
// 
//
//


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
///
/// note that this only includes mass and energy balances
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

/// for sonic flow, we need to get conditions where choked flow is 
/// achieved
pub mod choked_flow;
