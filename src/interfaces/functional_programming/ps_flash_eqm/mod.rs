/// checks if (p,s) point is within validity range
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
use uom::si::f64::*;

fn check_if_within_ps_validity_region(p: Pressure, s: SpecificHeatCapacity,){
    if is_outside_pressure_range(p) {
        panic!("p,s point is outside pressure range");
    };

    if is_below_isotherm_t_273_15(p, s) {
        panic!("p,s point below 273.15K");
    };
    if is_above_isotherm_t_1073_15(p, s) {
        panic!("p,s point above 1073.15K");
    };
}
