use uom::si::{f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

use crate::region_4_vap_liq_equilibrium::sat_pressure_4;

use super::pt_flash_eqm::FwdEqnRegion;

// allows the user to check which region one is in based on a ph flash
//
// note that ph flash does not work in region 5
pub(crate) fn ph_flash_region(p: Pressure, h: AvailableEnergy) -> FwdEqnRegion {

    check_if_within_ph_validity_region(p, h);

    // if inside validity range, then we will start partitioning
    // first, we check if pressure is smaller or greater than 16.529 MPa
    // this is saturation pressure at 623.15K 

    let t_for_pressure_boundary = ThermodynamicTemperature::new::<kelvin>(623.15);
    let p_boundary = sat_pressure_4(t_for_pressure_boundary);

    let is_p_below_16_529_mpa = p < p_boundary;

    // if p is below 16.529 mpa, then we use the eqns for below 16.529 mpa 

    if is_p_below_16_529_mpa {

        let is_region_1
            = is_ph_point_subcooled_liquid_region1_and_below_16_529_mpa(p, h);
        if is_region_1 {
            return FwdEqnRegion::Region1;
        };
        let is_region_2
            = is_ph_point_superheat_vap_region2_and_below_16_529_mpa(p, h);

        if is_region_2 {
            return FwdEqnRegion::Region2;
        };
        // else return region 4 
        return FwdEqnRegion::Region4;


    } 
    // if pressure is above 16.529 mpa,
    // then we have region 1,2,3 and 4
    // but we need to check the enthalpy as well because there are several 
    // regimes
    // this is the two phase region

    

    

    todo!()
}


/// panics if outside validity region
/// see page 38 top left
///
fn check_if_within_ph_validity_region(p: Pressure, h: AvailableEnergy,){
    if is_outside_pressure_range(p) {
        panic!("p,h point is outside pressure range");
    };

    if is_below_isotherm_t_273_15(p, h) {
        panic!("p,h point below 273.15K");
    };
    if is_above_isotherm_t_1073_15(p, h) {
        panic!("p,h point above 1073.15K");
    };
}

/// checks if (p,h) point is within validity range
/// of region 1 to 4
pub(crate) mod validity_range;
pub(crate) use validity_range::*;


/// this checks for boundary between single phase regions (1 to 3) 
/// and multiphase region 4
pub(crate) mod boundaries_from_single_phase_regions_to_region_4_multiphase;
pub(crate) use boundaries_from_single_phase_regions_to_region_4_multiphase::*;

/// this checks for boundary in between single phase regions (1 to 3) 
pub(crate) mod boundaries_between_single_phase_regions;
pub(crate) use boundaries_between_single_phase_regions::*;
