/// equations for isentropic nozzles, including choked flow 
/// at sonic speeds
pub mod isentropic_nozzles;
pub use isentropic_nozzles::*;


/// these contain equations for generator
/// where flux and stuff are used
#[allow(non_snake_case)]
pub mod generator;
pub use generator::*;

