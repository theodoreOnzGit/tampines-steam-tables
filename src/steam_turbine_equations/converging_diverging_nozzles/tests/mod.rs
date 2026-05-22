
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
