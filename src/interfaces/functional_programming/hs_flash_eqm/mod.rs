
use uom::si::f64::*;


use super::pt_flash_eqm::FwdEqnRegion;
// allows the user to check which region one is in based on a ph flash
//
// note that ph flash does not work in region 5
pub fn hs_flash_region(h: AvailableEnergy, s: SpecificHeatCapacity) -> FwdEqnRegion {

    todo!();
}

/// note:
/// (h,s) flashes along the isotherms 273.15K are not implemented 
/// for simplicity to avoid iterations
pub mod validity_range;
