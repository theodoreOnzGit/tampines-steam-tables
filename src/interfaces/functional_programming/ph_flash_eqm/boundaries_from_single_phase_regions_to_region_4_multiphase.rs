// boundary between region 1 (subcooled liq) and 4 (VLE) is 
// given by the saturated liq condition
//
// those are determined by the following sat temp and sat pressure 
// equations:


use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;

// now, we are given (p,h) points,
// this enthalpy h can be greater or less than the saturated liquid 
// enthalpy (region 1)
//
// and based on that, we can see if this is in region 1 

use uom::si::f64::*;
pub(crate) fn is_ph_point_subcooled_liquid_region1(p: Pressure, h: AvailableEnergy) -> bool {

    // first let's get saturated liquid enthalpy
    let sat_pressure_ref = p;
    let sat_temperature_ref = sat_temp_4(sat_pressure_ref);

    // now saturated liquid enthalpy given (p_sat,t_sat)
    let h_sat_liq = h_tp_1(sat_temperature_ref, sat_pressure_ref);

    // if enthalpy is less than saturated liquid enthalpy, then it is subcooled liq

    if h < h_sat_liq {
        return true;
    };

    return false;

}


