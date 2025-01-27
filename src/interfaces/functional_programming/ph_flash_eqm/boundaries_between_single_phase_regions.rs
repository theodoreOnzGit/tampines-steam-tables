use crate::constants::T_C_KELVIN;
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;
use crate::region_3_single_phase_plus_supercritical_steam::t_boundary_2_3;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;
use uom::si::f64::*;
use uom::si::pressure::megapascal;
use uom::si::thermodynamic_temperature::kelvin;

/// see page 39
pub(crate) fn is_ph_point_region_1_and_above_16_529_mpa(p: Pressure, h: AvailableEnergy) -> bool {

    // before anything, check if enthalpy is within enthalpy validity range 
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

    // now let's get the boundary line enthalpy 

    let h_boundary_line = h_tp_1(ref_temperature, p);

    // if enthalpy is more than this boundary, then it's not region 1 
    // i use greater than and NOT greater or equal to 
    // because the boundary line belongs to reg 1

    if h > h_boundary_line {
        return false;
    };

    // else it's region 1
    //

    return true;


}

/// also see page 39 
pub(crate) fn is_ph_point_region_2_and_above_16_529_mpa(p: Pressure,
    h: AvailableEnergy) -> bool {

    // before anything, check if enthalpy is within enthalpy validity range 
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

    // now let's get the boundary line enthalpy 
    // first, get the appropriate temperature Tb23 
    let t_boundary_b23 = t_boundary_2_3(p);

    let h_boundary_line = h_tp_2(t_boundary_b23, p);

    // if enthalpy is below this boundary line
    //
    // then its outside region 2
    //
    // i use smaller than and NOT smaller or equal to 
    // because the boundary line belongs to reg 2
    if h < h_boundary_line {
        return false;
    };

    // otherwise it's inside region 2
    return true;


}


/// also see page 39 
pub(crate) fn is_ph_point_region_3_and_above_critical_point(p: Pressure,
    h: AvailableEnergy) -> bool {

    // before anything, check if enthalpy is within enthalpy validity range 
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

    // now let's get the boundary line enthalpy 
    // first, get the appropriate temperature Tb23 
    // this is the upper bound
    let t_boundary_b23 = t_boundary_2_3(p);
    let h_boundary_line_23 = h_tp_2(t_boundary_b23, p);
    let t_boundary_isotherm = ThermodynamicTemperature::new::<kelvin>(623.15);
    let h_boundary_isotherm = h_tp_1(t_boundary_isotherm, p);

    // if enthalpy is outside this boundary line
    //
    // then its outside region 3
    //
    // i use smaller than and NOT smaller or equal to 
    // because the boundary line belongs to reg 2
    if h >= h_boundary_line_23 {
        return false;
    };

    if h <= h_boundary_isotherm {
        return false;
    };

    // otherwise it's inside region 3
    return true;


}


