
/// for subsonic flow, outlet pressure is higher than the throat pressure 
/// and the entire process is considered isentropic
/// 
pub mod cd_nozzle_subsonic_flow;


/// for perfectly expanded flow 
/// 
/// in this case, the flow is choked, so there is an isentropic 
/// process from 
pub mod diverging_nozzle_perfectly_expanded_supersonic;

/// for overexpanded flow, 
pub mod cd_nozzle_choked_flow_overexpanded;


/// Marviken tests 
pub mod marviken_tests;

/// From Figure 1 of:
///
/// Moody, F. J. (1975). Maximum discharge rate of liquid-vapor mixtures 
/// from vessels (No.
/// NEDO--21052). General Electric Co., San Jose, CA (United States). 
/// BWR Projects Dept..0 
///
/// Downloaded at:
/// https://www.osti.gov/servlets/purl/7309475
///
///
pub mod moody_critical_mass_flux_homogeneous_eqm;
