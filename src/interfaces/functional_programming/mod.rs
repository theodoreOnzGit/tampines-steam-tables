/// allows for pressure and temperature flash 
/// for all other properties 
/// (except steam quality, which cannot be 
/// determined via pt flashing)
///
/// this uses the forward equations
///
/// water/steam is assumed at thermodynamic equilibrium,
/// ie not metastable
pub mod pt_flash;

/// this pt_flash allows for metastable steam
/// which is NOT at thermodynamic equilibrium
pub mod pt_flash_metastable;
