// boundary between region 1 (subcooled liq) and 4 (VLE) is 
// given by the saturated liq condition
//
// those are determined by the following sat temp and sat pressure 
// equations:


use crate::backward_eqn_ph_region_1_to_4::p_s3_h;
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::kelvin;

/// see page 38 
pub(crate) fn is_ph_point_region_4_and_above_16_529_mpa(p: Pressure, h: AvailableEnergy) -> bool {

    // before anything, check if enthalpy is within enthalpy validity range 
    let ref_temperature = ThermodynamicTemperature::new::<kelvin>(623.15);
    let ref_pressure = sat_pressure_4(ref_temperature);

    let h_min_sat_liq = h_tp_1(ref_temperature, ref_pressure);
    let h_max_sat_vap = h_tp_2(ref_temperature, ref_pressure);

    if h > h_max_sat_vap {
        panic!(" enthalpy of p,h point is outside validity range");
    };
    if h < h_min_sat_liq {
        panic!(" enthalpy of p,h point is outside validity range");
    };

    // now, if within this range, we can check pressure
    // in comparison to the saturation enthalpy
    let mut p_sat_line = p_s3_h(h);



    // if the h is very close to the two phase region 
    // we need to correct this slightly
    //
    // This is in the section after table 2.30
    if ((p - p_sat_line)/p).get::<ratio>() < 1e-4 {
        p_sat_line = p_sat_line * (1.0 - 4.3e-6);
    };
    

    // if pressure is greater than this pressure, then it is region 3 
    // otherwise it's region 4 

    if p < p_sat_line {

        return true;

    };

    return false;

}
/// see page 38 
pub(crate) fn is_ph_point_region_3_and_from_16_529_mpa_to_crit_temp(p: Pressure, h: AvailableEnergy) -> bool {

    // before anything, check if enthalpy is within enthalpy validity range 
    let ref_temperature = ThermodynamicTemperature::new::<kelvin>(623.15);

    let t_low_bound = ThermodynamicTemperature::new::<kelvin>(623.15);

    let h_min = h_tp_1(t_low_bound, p);
    let h_max = h_tp_2(ref_temperature, p);

    if h > h_max {
        panic!(" enthalpy of p,h point is outside validity range");
    };
    if h < h_min {
        panic!(" enthalpy of p,h point is outside validity range");
    };

    todo!();

    // now, if within this range, we can check pressure
    // in comparison to the saturation enthalpy
    let mut p_sat_line = p_s3_h(h);



    // if the h is very close to the two phase region 
    // we need to correct this slightly
    //
    // This is in the section after table 2.30
    if ((p - p_sat_line)/p).get::<ratio>() < 1e-4 {
        p_sat_line = p_sat_line * (1.0 - 4.3e-6);
    };
    

    // if pressure is greater than this pressure, then it is region 3 
    // otherwise it's region 4 

    if p < p_sat_line {

        return true;

    };

    return false;

}


/// now, we are given (p,h) points,
/// this enthalpy h can be greater or less than the saturated liquid 
/// enthalpy (region 1)
///
/// and based on that, we can see if this is in region 1 
pub(crate) fn is_ph_point_subcooled_liquid_region1_and_below_16_529_mpa(p: Pressure, h: AvailableEnergy) -> bool {

    // before anything, check if it is within pressure validity range 
    let max_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(623.15));
    let min_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));

    if p > max_pressure || p < min_pressure {
        panic!(" pressure of p,h point is outside validity range");
    };

    // first let's get saturated liquid enthalpy
    let sat_temperature_ref = sat_temp_4(p);


    // now saturated liquid enthalpy given (p_sat,t_sat)
    let h_sat_liq = h_tp_1(sat_temperature_ref, p);

    // if enthalpy is less than saturated liquid enthalpy, then it is subcooled liq

    if h <= h_sat_liq {
        return true;
    };

    return false;

}

pub(crate) fn is_ph_point_superheat_vap_region2_and_below_16_529_mpa(p: Pressure, h: AvailableEnergy) -> bool {

    // before anything, check if it is within pressure validity range 
    let max_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(623.15));
    let min_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));

    if p > max_pressure || p < min_pressure {
        panic!(" pressure of p,h point is outside validity range");
    };
    // first let's get saturated vap enthalpy
    let sat_pressure_ref = p;
    let sat_temperature_ref = sat_temp_4(sat_pressure_ref);

    // now saturated liquid enthalpy given (p_sat,t_sat)
    let h_sat_vap = h_tp_2(sat_temperature_ref, sat_pressure_ref);

    // if enthalpy is less than saturated liquid enthalpy, then it is subcooled liq

    if h >= h_sat_vap {
        return true;
    };

    return false;

}



