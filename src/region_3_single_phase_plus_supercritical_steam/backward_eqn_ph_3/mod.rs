/// obtains temperature from pressure and enthalpy in region 3
pub mod t_ph_flash;
pub use t_ph_flash::*;

/// obtains specific volume from pressure and enthalpy in region 3
pub mod v_ph_flash;
pub use v_ph_flash::*;

/// boundary equations for 3a and 3b
pub mod boundary_eqn_3a_3b;
pub use boundary_eqn_3a_3b::*;
