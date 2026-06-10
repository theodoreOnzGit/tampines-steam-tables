
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

