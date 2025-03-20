
use uom::si::{f64::*,  specific_heat_capacity::kilojoule_per_kilogram_kelvin};
use validity_range::s_crit;


use super::pt_flash_eqm::FwdEqnRegion;
/// allows the user to check which region one is in based on a ph flash
///
/// note that ph flash does not work in region 5
///
/// the way to do region separation is first by entropy according to 
/// fig 2.14
///
/// once that is done, then we separate region by enthalpy.
pub fn hs_flash_region(h: AvailableEnergy, s: SpecificHeatCapacity) -> FwdEqnRegion {

    // this is absolute minimum and max entropy
    // based on fig 2.14
    let s_min = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(-0.00858);
    let s_max = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(11.92);

    if s < s_min {
        panic!("entropy too low for hs flash");
    };
    if s > s_max {
        panic!("entropy too high for hs flash");
    };

    let s_bound_low_entropy_region_1_and_4 = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.398);

    if s < s_bound_low_entropy_region_1_and_4 {
        hs_region_low_entropy_region_1_and_4(h, s);
    };

    let s_bound_low_entropy_region_1_3a_and_4 = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.778);

    if s <= s_bound_low_entropy_region_1_3a_and_4 {
        hs_region_low_entropy_region_1_3a_and_4(h, s);
    };

    let s_crit = s_crit();

    if s < s_crit {
        hs_region_near_crit_entropy_region_3a_and_4(h, s);
    };

    if s == s_crit {
        hs_region_crit_entropy_region_3a_and_4(h, s);
    };

    // this is for the region where TB23 and p2c are used
    let s_b23_min = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.048_096_828
        );

    let s_b23_max = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.260_578_707
        );

    if s < s_b23_min {
        hs_region_near_crit_entropy_region_3b_and_4(h, s);
    };

    if s <= s_b23_max {
        hs_region_near_crit_entropy_region_2c_3b_and_4(h, s);
    };
    
    let s_2bc = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.85
        );

    if s <= s_2bc {
        hs_region_high_entropy_region_2c_and_4(h, s);
    };

    // based on fig 2.19
    let s_where_2b_starts = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            6.070
        );

    if s < s_where_2b_starts {
        hs_region_high_entropy_region_2b_and_4(h, s);
    };

    // based on fig 2.19
    let s_where_2a_starts = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            7.852
        );

    if s <= s_where_2a_starts {
        hs_region_high_entropy_region_2b_2a_and_4(h, s);
    };

    // otherwise, should be in region 2a 
    
    hs_region_high_entropy_region_2a_only(h, s)

}

fn hs_region_low_entropy_region_1_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}

fn hs_region_low_entropy_region_1_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}

fn hs_region_near_crit_entropy_region_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_crit_entropy_region_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_near_crit_entropy_region_3b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_near_crit_entropy_region_2c_3b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_high_entropy_region_2c_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_high_entropy_region_2b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_high_entropy_region_2b_2a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}
fn hs_region_high_entropy_region_2a_only(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> FwdEqnRegion {

    todo!();

}

/// note:
/// (h,s) flashes along the isotherms 273.15K are not implemented 
/// for simplicity to avoid iterations
pub mod validity_range;
