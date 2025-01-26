#[warn(missing_docs)]
/// constants for the steam table calculations
pub mod constants;

/// region 1 
///
/// Temperature from 273.15 to 623.15 K 
/// Pressure from 0 to 100 MPa
///
/// Up to the saturation line. 
/// This I believe is subcooled liquid region
pub mod region_1_subcooled_liquid;



/// region 2 
///
/// vapour region
pub mod region_2_vapour;

/// region 3 
///
/// vapour liquid equilibrium
/// wet steam region, also includes supercritical region 
/// and critical point
///
/// auxilliary equation for region 2 and 3 are also put here
///
///
pub mod region_3_single_phase_plus_supercritical_steam;

/// region 4
///
/// two phase region
pub mod region_4_vap_liq_equilibrium;


/// region 5 
///
/// superheated steam region (ultra high temp)
pub mod region_5_superheated_steam;


/// backward equations ph boundary equations
/// overall equation
pub mod backward_eqn_ph_region_1_to_4;


/// public facing interfaces where the user 
/// simply inputs pressure and temperature 
/// or pressure and enthalpy etc 
/// and gets all the required data automatically
///
/// the logic for splitting between regions is 
/// mostly here 
pub mod interfaces;

/// allows for easy importing as with most rust 
/// crates. 
pub mod prelude;

