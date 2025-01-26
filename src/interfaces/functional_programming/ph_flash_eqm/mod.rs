use uom::si::{f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

use crate::region_4_vap_liq_equilibrium::sat_pressure_4;

use super::pt_flash_eqm::{h_tp_eqm, FwdEqnRegion};

// allows the user to check which region one is in based on a ph flash
//
// note that ph flash does not work in region 5
pub(crate) fn ph_flash_region(p: Pressure, h: AvailableEnergy) -> FwdEqnRegion {

    todo!()



}

// checks if pressure is 
// lower than saturation pressure at 273.15K or higher than 100 MPa
//
// if so, it falls outside the ph boundary
fn is_outside_pressure_range(p: Pressure) -> bool {

    // first determine if p,h point is outside pressure range
    let lower_pressure_limit: Pressure = 
        sat_pressure_4(
            ThermodynamicTemperature::new::<kelvin>(273.15)
        );

    let upper_pressure_limit: Pressure = Pressure::new::<megapascal>(100.0);

    if p < lower_pressure_limit {
        return true;
    };

    if p > upper_pressure_limit {
        return true;
    };

    return false;

}

// making a function to check if a p,h value is below the isotherm at 
// 273.15K
//
fn is_below_isotherm_t_273_15(p: Pressure, h: AvailableEnergy) -> bool{

    // first check if outside pressure range 
    if is_outside_pressure_range(p) {
        panic!("outside pressure range");
    };

    // let's have the lower enthalpy range 
    let lower_temp_bound = ThermodynamicTemperature::new::<kelvin>(273.15);

    let lower_bound_enthalpy = h_tp_eqm(lower_temp_bound, p);

    if h < lower_bound_enthalpy {
        return true;
    };

    return false;

}


