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
use uom::si::{f64::*, pressure::megapascal, ratio::ratio, thermodynamic_temperature::kelvin};

use crate::{region_1_subcooled_liquid::{s_tp_1, t_ps_1, u_tp_1, v_tp_1}, region_2_vapour::{s_tp_2, t_ps_2, u_tp_2, v_tp_2}, region_3_single_phase_plus_supercritical_steam::{s_3a3b_backwards_ps_boundary, s_rho_t_3, t_ps_3, u_rho_t_3, v_ps_3, v_tp_3c, v_tp_3r, v_tp_3s, v_tp_3t, v_tp_3u, v_tp_3x, v_tp_3y, v_tp_3z}, region_4_vap_liq_equilibrium::{sat_pressure_4, sat_temp_4}, region_5_steam_at_800_plus_degc::u_tp_5};

use super::pt_flash_eqm::FwdEqnRegion;
/// obtains temperature given pressure and entropy
pub fn t_ps_eqm(p: Pressure, s: SpecificHeatCapacity,) -> ThermodynamicTemperature {
    let region = ps_flash_region(p, s);

    match region {
        FwdEqnRegion::Region1 => t_ps_1(p, s),
        FwdEqnRegion::Region2 => t_ps_2(p, s),
        FwdEqnRegion::Region3 => t_ps_3(p, s),
        FwdEqnRegion::Region4 => {
            // if region 4, then just use the pressure to 
            // determine sat liq/vap temperature 
            sat_temp_4(p)
        },
        FwdEqnRegion::Region5 => todo!("region 5 ps flash not implemented"),
    }
}

/// obtains volume given pressure and entropy (except for region 5)
pub fn v_ps_eqm(p: Pressure, s: SpecificHeatCapacity) -> SpecificVolume {
    let region = ps_flash_region(p, s);

    match region {
        FwdEqnRegion::Region1 => {
            let t = t_ps_1(p, s);
            v_tp_1(t, p)
        },
        FwdEqnRegion::Region2 => {
            let t = t_ps_2(p, s);
            v_tp_2(t, p)
        },
        FwdEqnRegion::Region3 => v_ps_3(p, s),
        FwdEqnRegion::Region4 => {
            // in region 4 we get steam quality first
            // and then sat temp
            let steam_quality = x_ps_flash(p, s);
            let t_sat = sat_temp_4(p);

            // there are two cases for this... in region 4 
            // boundary one is where temperature is below 
            // 623.15 K 
            // and one where it is above 623.15 K 
            // in the first case, using region 1 and 2 is ok 
            // but in the latter case, it is not
            // one has to use region 3 eqns
            let t_sat_kelvin = t_sat.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let v_liq = v_tp_1(t_sat, p);
                let v_vap = v_tp_2(t_sat, p);

                let v = steam_quality * v_vap + (1.0 - steam_quality) * v_liq;

                v
            } else {
                // in the case we hit region 3 and region 4 boundary, the algorithm must differ

                // looks like all I need to do is use the backward equations 
                // v(t,p)
                // let's do the liquid one first 

                // this is with reference to the pt equations on 
                // fig 2.24 page 109
                //
                // alternatively, this is a cheat way... 
                // I look at t_sat, and force it to find liquid volume 
                // using a slightly colder temperature than tsat 
                // and for vapour I get it slightly higher than  
                // the tsat
                let p_mpa = p.get::<megapascal>();
                let v_vap: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 640.691 {
                        v_tp_3t(t_sat, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t_sat, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t_sat, p)
                    } else {
                        v_tp_3z(t_sat, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t_sat, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t_sat, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t_sat, p)
                    } else {
                        v_tp_3y(t_sat, p)
                    }


                };
                // yea this didn't work well... 
                // I'm going to use the tp equations in their individual regions

                let v = steam_quality * v_vap + (1.0 - steam_quality) * v_liq;

                v


            }
        },
        FwdEqnRegion::Region5 => todo!("ps flash not implemented for region 5"),
    }
}
/// obtains steam quality (vap fraction) given 
/// pressure and entropy 
pub fn x_ps_flash(p: Pressure, s: SpecificHeatCapacity,) -> f64 {

    let region = ps_flash_region(p, s);
    match region {
        // region 1 is liquid (but above crit point, doesn't really mater
        FwdEqnRegion::Region1 => 0.0,
        // region 2 is vapour, but above crit point doesn't really matter
        FwdEqnRegion::Region2 => 1.0,
        FwdEqnRegion::Region3 => {
            // region 3 is special, if it is equal or above 
            // crit point, then just consider it vapour,
            // doesn't really matter
            if p >= Pressure::new::<megapascal>(22.064) {
                return 1.0;
            };

            // otherwise, in region 3, we check if 
            // the entropy is below critical entropy 
            let s_crit = s_3a3b_backwards_ps_boundary();

            if s < s_crit {
                return 0.0;
            } else {
                return 1.0;
            };
        },
        FwdEqnRegion::Region4 => {
            // for this we consider vapour liquid equilibrium
            //
            // h = hvap (x) + hliq (1-x)
            // h = x (hvap - hliq) + hliq
            // h-hliq = x (hvap - hliq)
            // x = (h-hliq)/(hvap - hliq)
            
            // first let's take saturation pressure 
            let t_sat = sat_temp_4(p);


            let t_sat_kelvin = t_sat.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let s_liq = s_tp_1(t_sat, p);
                let s_vap = s_tp_2(t_sat, p);

                let x: Ratio = (s-s_liq)/(s_vap-s_liq);
                return x.get::<ratio>();
            } else {
                // in the case we hit region 3 and region 4 boundary, the algorithm must differ

                // looks like all I need to do is use the backward equations 
                // v(t,p)
                // let's do the liquid one first 

                // this is with reference to the pt equations on 
                // fig 2.24 page 109
                //
                // alternatively, this is a cheat way... 
                // I look at t_sat, and force it to find liquid volume 
                // using a slightly colder temperature than tsat 
                // and for vapour I get it slightly higher than  
                // the tsat
                let v_vap: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 640.691 {
                        v_tp_3t(t_sat, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t_sat, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if t_sat_kelvin <= 646.599 {
                        v_tp_3x(t_sat, p)
                    } else {
                        v_tp_3z(t_sat, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t_sat, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t_sat, p)
                    } else if t_sat_kelvin <= 646.483 {
                        v_tp_3u(t_sat, p)
                    } else {
                        v_tp_3y(t_sat, p)
                    }


                };
                // now we have t v equations, we can get h_liq 
                // and h_vap
                // using t,v equations

                let h_liq = s_rho_t_3(v_liq.recip(), t_sat);
                let h_vap = s_rho_t_3(v_vap.recip(), t_sat);

                if h_liq == h_vap {
                    // at supercritical point, just assume vapour
                    // this prevents a divide by zero error
                    return 1.0;
                };


                let x: Ratio = (s-h_liq)/(h_vap-h_liq);
                return x.get::<ratio>();


            }

        },
        // this is a placeholder, 
        // but technically if in region 5, the vapour 
        // quality is 1.0
        FwdEqnRegion::Region5 => 1.0,
    }

}


/// returns the internal energy given entropy and pressure
pub fn u_ph_eqm(p: Pressure, s: SpecificHeatCapacity) -> AvailableEnergy {
    let t = t_ps_eqm(p, s);
    let region = ps_flash_region(p, s);

    match region {
        FwdEqnRegion::Region1 => u_tp_1(t, p),
        FwdEqnRegion::Region2 => u_tp_2(t, p),
        FwdEqnRegion::Region3 => {
            // near supercritical point, we have to be a little bit more careful
            // need to take into account steam quality too

            // so first, test if we are JUST on the saturation line
            // because that is part of the check
            //
            // or rather, just get the volume first
            // this will make it much easier
            let v = v_ps_eqm(p, s);
            let rho = v.recip();
            u_rho_t_3(rho, t)
        },
        FwdEqnRegion::Region4 => {
            // for region 4 specifically, we determine 
            // steam quality first
            let t_sat = sat_temp_4(p);
            let t_sat_kelvin = t_sat.get::<kelvin>();
            let steam_quality = x_ps_flash(p, s);
            if t_sat_kelvin <= 623.15 {

                let u_liq = u_tp_1(t_sat, p);
                let u_vap = u_tp_2(t_sat, p);

                let u = steam_quality * u_vap + (1.0 - steam_quality) * u_liq;

                u
            } else {
                // in the case we hit region 3, the algorithm must differ

                // looks like all I need to do is use the backward equations 
                // v(t,p)
                // let's do the liquid one first 

                // this is with reference to the pt equations on 
                // fig 2.24 page 109
                //
                // alternatively, this is a cheat way... 
                // I look at t_sat, and force it to find liquid volume 
                // using a slightly colder temperature than tsat 
                // and for vapour I get it slightly higher than  
                // the tsat

                let p_mpa = p.get::<megapascal>();
                let v_vap: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 640.691 {
                        v_tp_3t(t_sat, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t_sat, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t_sat, p)
                    } else {
                        v_tp_3z(t_sat, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t_sat, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t_sat, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t_sat, p)
                    } else {
                        v_tp_3y(t_sat, p)
                    }


                };

                // now let's get u for liquid and vapour 

                let u_liq = u_rho_t_3(v_liq.recip(), t);
                let u_vap = u_rho_t_3(v_vap.recip(), t);


                let u = steam_quality * u_vap + (1.0 - steam_quality) * u_liq;

                u

            }
        },
        FwdEqnRegion::Region5 => u_tp_5(t, p),
    }
}

// allows the user to check which region one is in based on a ps flash
//
// note that ps flash does not work in region 5
pub fn ps_flash_region(p: Pressure, s: SpecificHeatCapacity) -> FwdEqnRegion {

    check_if_within_ps_validity_region(p, s);

    // if inside validity range, then we will start partitioning
    // first, we check if pressure is smaller or greater than 16.529 MPa
    // this is saturation pressure at 623.15K 

    let t_for_pressure_boundary = ThermodynamicTemperature::new::<kelvin>(623.15);
    let p_boundary = sat_pressure_4(t_for_pressure_boundary);

    let is_p_below_16_529_mpa = p < p_boundary;

    // if p is below 16.529 mpa, then we use the eqns for below 16.529 mpa 

    if is_p_below_16_529_mpa {

        let is_region_1
            = is_ps_point_subcooled_liquid_region1_and_below_16_529_mpa(p, s);
        if is_region_1 {
            return FwdEqnRegion::Region1;
        };
        let is_region_2
            = is_ps_point_superheat_vap_region2_and_below_16_529_mpa(p, s);

        if is_region_2 {
            return FwdEqnRegion::Region2;
        };
        // else return region 4 
        return FwdEqnRegion::Region4;


    } 
    // if pressure is above 16.529 mpa,
    // then we have region 1,2,3 and 4
    // but we need to check the entropy as well because there are several 
    // regimes this is the two psase region
    // this is from h = 1670.9 kJ/kg to about 2563.6 kJ/kg
    // and from 16.529 MPa to 22.064 Mpa (critical point)
    // here is where things are potentially two psase

    let is_region_1 
        = is_ps_point_region_1_and_above_16_529_mpa(p, s);

    if is_region_1 {
        return FwdEqnRegion::Region1;
    };

    let is_region_2 
        = is_ps_point_region_2_and_above_16_529_mpa(p, s);
    if is_region_2 {
        return FwdEqnRegion::Region2;
    };

    // the checks for if it is region 1 or 2 already exclude points 
    // outside h = 1670.9 kJ/kg to about 2563.6 kJ/kg
    //
    // this is above 22.064 MPa

    let is_region_3_above_supercrit 
        = is_ps_point_region_3_and_above_critical_point(p, s);

    if is_region_3_above_supercrit {
        return FwdEqnRegion::Region3;
    };

    // we want to check if this is region 3 and in the range 
    // 16.520 MPa up to crit point 22.064 MPa
    let is_region_3_below_supercrit
        = is_ps_point_region_3_and_from_16_529_mpa_to_crit_temp(p, s);

    if is_region_3_below_supercrit {
        return FwdEqnRegion::Region3;
    };

    // now we shall have to decide if it is region 3 or 4

    let is_region_4 = 
        is_ps_point_region_4_and_above_16_529_mpa(p, s);
    if is_region_4 {
        return FwdEqnRegion::Region4;
    };

    

    // otherwise it's region 3 

    return FwdEqnRegion::Region3;
}
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
