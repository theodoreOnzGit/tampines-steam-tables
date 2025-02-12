use crate::constants::T_C_KELVIN;
use crate::region_1_subcooled_liquid::s_tp_1;
use crate::region_2_vapour::s_tp_2;
use crate::region_3_single_phase_plus_supercritical_steam::t_boundary_2_3;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;
use uom::si::f64::*;
use uom::si::pressure::megapascal;
use uom::si::thermodynamic_temperature::kelvin;

/// see page 55-56
pub(crate) fn is_ps_point_region_1_and_above_16_529_mpa(
    p: Pressure, s: SpecificHeatCapacity) -> bool {

    // before anything, check if entropy is within entropy validity range 
    let ref_temperature = ThermodynamicTemperature::new::<kelvin>(623.15);
    let min_pressure = sat_pressure_4(ref_temperature);
    let max_pressure = Pressure::new::<megapascal>(100.0);

    if p < min_pressure {
        panic!("p in (p,h) point is outside validity range");
    };
    if p > max_pressure {
        panic!("p in (p,h) point is outside validity range");
    };

    // note that points along this boundary belongs to region 1 
    // (see page 11)
    // also, points along the P_B23 boundary line belong to region 2

    // now let's get the boundary line entropy 

    let s_boundary_line = s_tp_1(ref_temperature, p);

    // if entropy is more than this boundary, then it's not region 1 
    // i use greater than and NOT greater or equal to 
    // because the boundary line belongs to reg 1

    if s > s_boundary_line {
        return false;
    };

    // else it's region 1
    //

    return true;


}

/// also see page 55-56
pub(crate) fn is_ps_point_region_2_and_above_16_529_mpa(p: Pressure,
    s: SpecificHeatCapacity) -> bool {

    // before anything, check if entropy is within entropy validity range 
    let ref_temperature = ThermodynamicTemperature::new::<kelvin>(623.15);
    let min_pressure = sat_pressure_4(ref_temperature);
    let max_pressure = Pressure::new::<megapascal>(100.0);

    if p < min_pressure {
        panic!("p in (p,h) point is outside validity range");
    };
    if p > max_pressure {
        panic!("p in (p,h) point is outside validity range");
    };

    // note that points along this boundary belongs to region 1 
    // (see page 11)
    // also, points along the P_B23 boundary line belong to region 2

    // now let's get the boundary line entropy 
    // first, get the appropriate temperature Tb23 
    let t_boundary_b23 = t_boundary_2_3(p);

    let s_boundary_line = s_tp_2(t_boundary_b23, p);

    // if entropy is below this boundary line
    //
    // then its outside region 2
    //
    // i use smaller than and NOT smaller or equal to 
    // because the boundary line belongs to reg 2
    if s < s_boundary_line {
        return false;
    };

    // otherwise it's inside region 2
    return true;


}


/// also see page 39 
pub(crate) fn is_ps_point_region_3_and_above_critical_point(p: Pressure,
    s: SpecificHeatCapacity) -> bool {

    // before anything, check if entropy is within entropy validity range 
    let ref_temperature = ThermodynamicTemperature::new::<kelvin>(T_C_KELVIN);
    let min_pressure = sat_pressure_4(ref_temperature);
    let max_pressure = Pressure::new::<megapascal>(100.0);

    if p < min_pressure {
        return false;
    };
    if p > max_pressure {
        panic!("p in (p,h) point is outside validity range");
    };

    // note that points along this boundary belongs to region 1 
    // (see page 11)
    // also, points along the P_B23 boundary line belong to region 2

    // now let's get the boundary line entropy 
    // first, get the appropriate temperature Tb23 
    // this is the upper bound
    let t_boundary_b23 = t_boundary_2_3(p);
    let s_boundary_line_23 = s_tp_2(t_boundary_b23, p);
    let t_boundary_isotherm = ThermodynamicTemperature::new::<kelvin>(623.15);
    let s_boundary_isotherm = s_tp_1(t_boundary_isotherm, p);

    // if entropy is outside this boundary line
    //
    // then its outside region 3
    //
    // i use smaller than and NOT smaller or equal to 
    // because the boundary line belongs to reg 2
    if s >= s_boundary_line_23 {
        return false;
    };

    if s <= s_boundary_isotherm {
        return false;
    };

    // otherwise it's inside region 3
    return true;


}



