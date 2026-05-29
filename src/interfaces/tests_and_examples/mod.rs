
/// these are tests to check the functionality 
/// of ph flash regions
pub mod ph_flash_regions;

/// these are tests to check the functionality 
/// of hs flash regions
/// note: does not include out of bounds just yet..
pub mod hs_flash_regions;

/// aims to reproduce steam tables using ph flash
#[cfg(test)]
pub mod ph_flash_steam_table;

/// aims to reproduce steam tables using pt flash
#[cfg(test)]
pub mod pt_flash_steam_table;
/// aims to reproduce steam tables using ps flash
#[cfg(test)]
pub mod ps_flash_steam_table;
///// aims to reproduce steam tables using hs flash
#[cfg(test)]
pub mod hs_flash_steam_table;

/// this aims to produce critical pressure values for 
/// VLE 
/// Based on figure 2 of:
/// https://www.osti.gov/servlets/purl/7309475
///
/// for speed of sound in VLE mixtures, refer to:
///
/// https://geology.illinois.edu/~skieffer/papers/SoundSpeed_JGR1977.pdf
#[cfg(test)]
pub mod critical_pressure_moody_fig2;

/// Saha, P. (1978). A review of two-phase steam-water 
/// critical flow models with emphasis on thermal nonequilibrium.
/// https://www.nrc.gov/docs/ML1925/ML19256F779.pdf
///
/// This provides a homogeneous equilibrium model (HEM) for critical 
/// flow
///
/// From page 2-5 of Saha's publication.
///
/// G = rho_mean * u_mean 
///
/// h_0 = h + 0.5 * u_mean * u_mean 
///
/// For homogeneous flow, 
/// 1/rho_m = (1-x)/rho_l + x/rho_v
///
/// h = (1-x) h_l + x h_v 
///
/// Assuming flow comes out saturated at the critical Pressure (P):
///
/// T_l = T_v = T_sat (P)
///
/// The vapour and liquid properties then take on their saturated 
/// properties
/// rho_l = rho_l (P) 
/// rho_v = rho_v (P)
/// h_l = h_l (P) 
/// h_v = h_v (P)
///
/// We then substitute these values to find u_mean, and rho_mean
///
/// The quality x as a function of the critical pressure (a saturation 
/// pressure) is:
///
/// x = (s_0 - s_l(P))/(s_v (P) - s_l (P))
///
/// Of course, changing P will change x, and we are just changing the 
/// pressure values to find the maximum flowrate. This is an optimisation 
/// problem.
///
/// For this, the critical pressure and critical mass flux are returned 
/// as a pair. Nothing here deals with sonic velocity. So speed of sound 
/// is not really a matter here.
///
/// There also needs to be an algorithm with which to obtain these
/// minimum points outside the speed of sound.
/// 
///
#[cfg(test)]
pub mod critical_flow_hem;

