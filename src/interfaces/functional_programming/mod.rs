/// allows for pressure and temperature flash 
/// for all other properties 
/// (except steam quality, which cannot be 
/// determined via pt flashing)
///
/// this uses the forward equations
///
/// water/steam is assumed at thermodynamic equilibrium,
/// ie not metastable
pub mod pt_flash_eqm;

/// allows for pressure enthalpy flash
pub mod ph_flash_eqm;

/// this pt_flash allows for metastable steam
/// which is NOT at thermodynamic equilibrium
///
/// this mostly deals with areas around region 2
pub mod pt_flash_metastable;


/// allows for pressure entropy flash 
pub mod ps_flash_eqm;
