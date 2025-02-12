// boundary between region 1 (subcooled liq) and 4 (VLE) is 
// given by the saturated liq condition
//
// those are determined by the following sat temp and sat pressure 
// equations:


use crate::backward_eqn_ps_region_1_to_4::boundary_eqn_ps3::p_s3_s;
use crate::region_1_subcooled_liquid::s_tp_1;
use crate::region_2_vapour::s_tp_2;
use crate::region_3_single_phase_plus_supercritical_steam::t_boundary_2_3;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::kelvin;

/// see page 55
pub(crate) fn is_ph_point_region_4_and_above_16_529_mpa(
    p: Pressure, s: SpecificHeatCapacity) -> bool {

    // before anything, check if entropy is within entropy validity range 
    let ref_temperature = ThermodynamicTemperature::new::<kelvin>(623.15);
    let ref_pressure = sat_pressure_4(ref_temperature);

    let s_min_sat_liq = s_tp_1(ref_temperature, ref_pressure);
    let s_max_sat_vap = s_tp_2(ref_temperature, ref_pressure);

    if s > s_max_sat_vap {
        panic!(" entropy of p,s point is outside validity range");
    };
    if s < s_min_sat_liq {
        panic!(" entropy of p,s point is outside validity range");
    };

    // now, if within this range, we can check pressure
    // in comparison to the saturation entropy
    let mut p_sat_line = p_s3_s(s);



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
/// see page 55
pub(crate) fn is_ps_point_region_3_and_from_16_529_mpa_to_crit_temp(
    p: Pressure, s: SpecificHeatCapacity) -> bool {

    // before anything, check if entropy is within entropy validity range 
    let t_high_bound = t_boundary_2_3(p);
    let t_low_bound = ThermodynamicTemperature::new::<kelvin>(623.15);

    let s_min = s_tp_1(t_low_bound, p);
    let s_max = s_tp_2(t_high_bound, p);

    if s > s_max {
        panic!("entropy of p,s point is outside validity range");
    };
    if s < s_min {
        panic!("entropy of p,s point is outside validity range");
    };

    // next is to see the temperature

    // now, if within this range, we can check pressure
    // in comparison to the saturation entropy
    let mut p_sat_line = p_s3_s(s);

    // if the h is very close to the two phase region 
    // we need to correct this slightly
    //
    // This is in the section after table 2.30
    if ((p - p_sat_line)/p).get::<ratio>() < 1e-4 {
        p_sat_line = p_sat_line * (1.0 - 4.3e-6);
    };
    

    // if pressure is greater than this pressure, then it is region 3 
    // otherwise it's region 4 

    if p >= p_sat_line {

        return true;

    };

    return false;

}


/// now, we are given (p,s) points,
/// this entropy s can be greater or less than the saturated liquid 
/// entropy (region 1)
///
/// and based on that, we can see if this is in region 1 
pub(crate) fn is_ps_point_subcooled_liquid_region1_and_below_16_529_mpa(
    p: Pressure, s: SpecificHeatCapacity) -> bool {

    // before anything, check if it is within pressure validity range 
    let max_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(623.15));
    let min_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));

    if p > max_pressure || p < min_pressure {
        panic!(" pressure of p,s point is outside validity range");
    };

    // first let's get saturated liquid entropy
    let sat_temperature_ref = sat_temp_4(p);


    // now saturated liquid entropy given (p_sat,t_sat)
    let s_sat_liq = s_tp_1(sat_temperature_ref, p);

    // if entropy is less than saturated liquid entropy, then it is subcooled liq

    if s <= s_sat_liq {
        return true;
    };

    return false;

}

pub(crate) fn is_ps_point_superheat_vap_region2_and_below_16_529_mpa(
    p: Pressure, s: SpecificHeatCapacity) -> bool {

    // before anything, check if it is within pressure validity range 
    let max_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(623.15));
    let min_pressure = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));

    if p > max_pressure || p < min_pressure {
        panic!(" pressure of p,s point is outside validity range");
    };
    // first let's get saturated vap entropy
    let sat_pressure_ref = p;
    let sat_temperature_ref = sat_temp_4(sat_pressure_ref);

    // now saturated liquid entropy given (p_sat,t_sat)
    let s_sat_vap = s_tp_2(sat_temperature_ref, sat_pressure_ref);

    // if entropy is less than saturated liquid entropy, then it is subcooled liq

    if s >= s_sat_vap {
        return true;
    };

    return false;

}




