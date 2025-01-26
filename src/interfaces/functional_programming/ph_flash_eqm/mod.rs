use uom::si::f64::*;

use super::pt_flash_eqm::FwdEqnRegion;

// allows the user to check which region one is in based on a ph flash
//
// note that ph flash does not work in region 5
pub(crate) fn ph_flash_region(p: Pressure, h: AvailableEnergy) -> FwdEqnRegion {

    todo!()



}

// making a function to check if a p,h value is below the isotherm at 
// 273.15K
fn is_below_isotherm_t_273_15(p: Pressure, h: AvailableEnergy) -> bool{

    todo!()
}
