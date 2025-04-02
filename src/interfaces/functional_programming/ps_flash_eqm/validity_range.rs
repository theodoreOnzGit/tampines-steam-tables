use uom::si::{f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

use crate::{interfaces::functional_programming::pt_flash_eqm::s_tp_eqm_single_phase, region_4_vap_liq_equilibrium::sat_pressure_4};

// checks if pressure is 
// lower than saturation pressure at 273.15K or higher than 100 MPa
//
// if so, it falls outside the ps boundary
pub(crate) fn is_outside_pressure_range(p: Pressure) -> bool {

    // first determine if p,h point is outside pressure range
    let lower_pressure_limit: Pressure = 
        sat_pressure_4(
            ThermodynamicTemperature::new::<kelvin>(273.15)
        );

    let upper_pressure_limit: Pressure = Pressure::new::<megapascal>(100.0);

    if p < lower_pressure_limit {
        dbg!(&(p,lower_pressure_limit));
        panic!("p,s point is lower than acceptable pressure range");
    };

    if p > upper_pressure_limit {
        dbg!(&(p,upper_pressure_limit));
        panic!("p,s point is higher than acceptable pressure range");
    };

    return false;

}

// making a function to check if a p,s value is below the isotherm at 
// 273.15K
//
pub(crate) fn is_below_isotherm_t_273_15(p: Pressure, s: SpecificHeatCapacity) -> bool{

    // first check if outside pressure range 
    if is_outside_pressure_range(p) {
        panic!("outside pressure range");
    };

    // let's have the lower enthalpy range 
    let lower_temp_bound = ThermodynamicTemperature::new::<kelvin>(273.15);

    let lower_bound_entropy = s_tp_eqm_single_phase(lower_temp_bound, p);

    if s < lower_bound_entropy {
        return true;
    };

    return false;

}

// making a function to check if p,h value is above the isotherm T = 1073.15K
pub(crate) fn is_above_isotherm_t_1073_15(p: Pressure,s: SpecificHeatCapacity) -> bool {
    // first check if outside pressure range 
    if is_outside_pressure_range(p) {
        panic!("outside pressure range");
    };

    let upper_temp_bound = ThermodynamicTemperature::new::<kelvin>(1073.15);

    let upper_bound_entropy = s_tp_eqm_single_phase(upper_temp_bound, p);

    if s > upper_bound_entropy {
        dbg!(&(s,upper_bound_entropy));
        return true;
    };

    return false;

}




