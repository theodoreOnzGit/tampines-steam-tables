use uom::si::f64::*;
use uom::si::volume::cubic_meter;

use crate::prelude::{TampinesSteamTableCV};
use crate::steam_turbine_equations::diverging_nozzle::guess_velocity_and_state_for_diverge_nozzle_from_choked_throat;
use crate::steam_turbine_equations::isentropic_converging_nozzle::get_choked_flow_massrate_and_state_from_stagnation_properties_and_area;

/// here is the main function meant to calculate mass flowrates between 
/// two control volumes using a converging diverging nozzle 
/// This is useful for a stator (CD nozzle) and impulse turbine section 
/// where pressure drop across this impulse turbine is (ideally) zero
///
/// given an inlet pressure and enthalpy, and velocity
/// flow is accelerated isentropically to the choke point
///
/// The control volume pressure and enthalpy are reflected by 
/// state 1 and state 2 below
/// (p1, h1) -> nozzle -> (p2, h2)
///
/// the outlet pressure is known (equal to p2), but not the outlet enthalpy 
/// and velocity. This will be calculated here
/// 
///
/// The algorithm here obtains the stagnation properties of the inlet 
/// by isentropically increasing enthalpy. 
///
/// It then calculates if choked flow occurs.
///
/// In the case of subsonic flow:
/// Then if no choked flow occurs, mass flowrates are based on the 
/// energy balance depending on nozzle dimensions. Isentropy is also assumed
///
///
/// note that for this code to work, p1 needs to be greater than p2
///
/// note for this, velocity v1 is not used to calculate mass flowrate, 
/// but just to get stagnation enthalpy
///
/// For subsonic flows, isentropy is assumed 
/// For sonic flows without chokes, isentropy is also assumed 
///
/// For sonic flows with underexpansion, Joule Thomson throttling is assumed 
/// as a conservative estimate for entropy generation
///
/// For sonic flows with over expansion, a (p,h) algorithm is used to 
/// iteratively determine the outlet flow properties.
#[inline]
pub fn calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity,
    a_throat: Area,
    a2: Area,
    p2: Pressure,
    ) -> (Velocity, MassRate, TampinesSteamTableCV) {

    assert!(p1 > p2);

    // first let's obtain stagnation enthalpy 

    let h0: AvailableEnergy = h1 + 0.5 * v1 * v1;
    // we get the entropy 
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();

    // now stagnation properties come from s1 = s0 
    let s0 = s1;
    let state_0 = TampinesSteamTableCV::new_from_hs(h0, s0, ref_vol);
    let p0 = state_0.get_pressure();

    // now we can calculate choked flow 
    // this assumes that pressure does go to choked flow state
    let (choked_mass_flowrate, choked_state) = 
        get_choked_flow_massrate_and_state_from_stagnation_properties_and_area(
            p0, h0, a_throat
        );
    let p_throat_critical = choked_state.get_pressure();

    // now for subsonic flow, we can just take pressure difference between 
    // these two states

    let s2_ideal = s1;
    let outlet_flow_state_ideal = 
        TampinesSteamTableCV::new_from_ps(p2, s2_ideal, ref_vol);
    let rho_out_ideal = outlet_flow_state_ideal.get_rho();
    let h_out_ideal = outlet_flow_state_ideal.get_specific_enthalpy();
    let v_out_ideal: Velocity = (2.0 * (h0 - h_out_ideal)).sqrt();
    let m_ideal: MassRate = rho_out_ideal * v_out_ideal * a2;

    let throat_has_choked_flow: bool = p2 < p_throat_critical;

    if !throat_has_choked_flow {
        // in this case, we have subsonic flow
        // we use a (p,s) algorithm to find the flowrate at the back
        // these are already calculated, we can just return this

        let state_outlet = 
            outlet_flow_state_ideal;

        let m = m_ideal;
        let v = v_out_ideal;

        return (v, m , state_outlet);


    }

    // this is in the case we do have choked flow
    // mass flowrate is fixed 
    let m = choked_mass_flowrate;
    // outlet flow is now decided upon using a (p,h) algorithm
    let state_throat = choked_state;

    let (v, state_outlet) = 
        guess_velocity_and_state_for_diverge_nozzle_from_choked_throat(
            h0, 
            s0,
            p2, 
            a_throat, 
            a2,
            m, 
            state_throat
        );




    return (v, m , state_outlet);
}

#[cfg(test)]
mod tests;


/// these are converging nozzles,
/// for this, usually we assume isentropy.
/// This part is relatively straightforward
pub mod isentropic_converging_nozzle;

/// this is the diverging part of the C-D nozzle
/// where (p,h) flashing is used instead of (p,s) or (h,s) flashing
/// shocks are not explictly calculated, but entropy increase is assumed
pub mod diverging_nozzle;

/// this covers differential nozzle equations,
/// which may not be accurate for all sonic and transonic flows
/// (incomplete, not production ready)
mod differential_nozzle_equations;




/// momentum balance for Rayleigh line (incomplete)
mod momentum_balance_rayleigh_line;

/// for sonic flow, we need to get conditions where choked flow is 
/// achieved
///
/// these are for textbook questions, where basic verification is performed
/// to see if choked flow calculation is correct
pub mod choked_flow;
