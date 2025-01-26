// between region 1 (liquid) and 4 (Vap liq eqm), the boundary is 
// the saturated liquid line 
// at  T = 273.15 K (freezing pt) to 
// T = 623.15 K (critical pt)
//
// between region 3 (supercritical-ish) and 4 (vap liq eqm) 
// the isotherm at critical temp 623.15 K is used 
// the small range of enthalpies is given by the 
// enthalpy of region 1 at sat liq 
// and enthalpy of region 2 at sat vap 
//
// between region 2 (vap) and 4 (vap liq eqm),
// the sat vapour line (dew pt) is from 
// the sat vapour eqn for region 2 from 273.15 K to 
// 673.15 K
//
// for boundary between 1 and 3, 
//
// the isotherm at T = 623.15 K is used (crit temp) 
// 
// the boundary at region 1 to 3 is part of region 1,
// and the boundary at region 2 to 3 is considered 
// region 2
//
// 
// for boundary between 2 and 3, the T-B23 line is used 
// 

pub mod boundary_eqn_ps3;
pub use boundary_eqn_ps3::*;

#[cfg(test)]
mod tests;
