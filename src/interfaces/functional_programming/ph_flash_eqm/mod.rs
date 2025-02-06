use uom::si::{f64::*, pressure::megapascal, ratio::ratio, thermodynamic_temperature::kelvin};

use crate::region_5_steam_at_800_plus_degc::*;
use crate::region_4_vap_liq_equilibrium::*;
use crate::region_3_single_phase_plus_supercritical_steam::*;
use crate::region_2_vapour::*;
use crate::region_1_subcooled_liquid::*;
use crate::region_1_subcooled_liquid::InversePressure;

use super::pt_flash_eqm::FwdEqnRegion;

/// obtains temperature given pressure and enthalpy
pub fn t_ph_eqm(p: Pressure, h: AvailableEnergy,) -> ThermodynamicTemperature {
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => t_ph_1(p, h),
        FwdEqnRegion::Region2 => t_ph_2(p, h),
        FwdEqnRegion::Region3 => t_ph_3(p, h),
        FwdEqnRegion::Region4 => {
            // if region 4, then just use the pressure to 
            // determine sat liq/vap temperature 
            sat_temp_4(p)
        },
        FwdEqnRegion::Region5 => todo!("region 5 ph flash not implemented"),
    }
}

/// obtains volume given pressure and enthalpy (except for region 5)
pub fn v_ph_eqm(p: Pressure, h: AvailableEnergy) -> SpecificVolume {
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => {
            let t = t_ph_1(p, h);
            v_tp_1(t, p)
        },
        FwdEqnRegion::Region2 => {
            let t = t_ph_2(p, h);
            v_tp_2(t, p)
        },
        FwdEqnRegion::Region3 => v_ph_3(p, h),
        FwdEqnRegion::Region4 => {
            // in region 4 we get steam quality first
            // and then sat temp
            let steam_quality = x_ph_flash(p, h);
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
        FwdEqnRegion::Region5 => todo!("ph flash not implemented for region 5"),
    }
}
/// returns the internal energy given temperature and pressure
pub fn u_ph_eqm(p: Pressure, h: AvailableEnergy) -> AvailableEnergy {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

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
            let v = v_ph_eqm(p, h);
            let rho = v.recip();
            u_rho_t_3(rho, t)
        },
        FwdEqnRegion::Region4 => {
            // for region 4 specifically, we determine 
            // steam quality first
            let t_sat = sat_temp_4(p);
            let t_sat_kelvin = t_sat.get::<kelvin>();
            let steam_quality = x_ph_flash(p, h);
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


/// returns the specific entropy given temperature and pressure
pub fn s_ph_eqm(p: Pressure, h: AvailableEnergy) -> SpecificHeatCapacity {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => s_tp_1(t, p),
        FwdEqnRegion::Region2 => s_tp_2(t, p),
        FwdEqnRegion::Region3 => {
            // near supercritical point, we have to be a little bit more careful
            // need to take into account steam quality too

            // so first, test if we are JUST on the saturation line
            // because that is part of the check
            //
            // or rather, just get the volume first
            // this will make it much easier
            let v = v_ph_eqm(p, h);
            let rho = v.recip();

            s_rho_t_3(rho, t)
        },
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let t_sat = sat_temp_4(p);
            let t_sat_kelvin = t_sat.get::<kelvin>();
            let steam_quality = x_ph_flash(p, h);

            if t_sat_kelvin <= 623.15 {

                let s_liq = s_tp_1(t_sat, p);
                let s_vap = s_tp_2(t_sat, p);

                let s = steam_quality * s_vap + (1.0 - steam_quality) * s_liq;

                s 
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

                let s_liq = s_rho_t_3(v_liq.recip(), t);
                let s_vap = s_rho_t_3(v_vap.recip(), t);

                let s = steam_quality * s_vap + (1.0 - steam_quality) * s_liq;

                s 

            }
        },
        FwdEqnRegion::Region5 => s_tp_5(t, p),
    }
}

/// returns the isobaric (const pressure) heat capacitygiven temperature and pressure
pub fn cp_ph_eqm(p: Pressure, h: AvailableEnergy) -> SpecificHeatCapacity {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => cp_tp_1(t, p),
        FwdEqnRegion::Region2 => cp_tp_2(t, p),
        FwdEqnRegion::Region3 => {

            // near supercritical point, we have to be a little bit more careful
            // need to take into account steam quality too

            // so first, test if we are JUST on the saturation line
            // because that is part of the check
            //
            // or rather, just get the volume first
            // this will make it much easier
            let v = v_ph_eqm(p, h);
            let rho = v.recip();
            cp_rho_t_3(rho, t)
        },
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let cp_liq = cp_tp_1(t_sat, p);
            let cp_vap = cp_tp_2(t_sat, p);

            let cp = steam_quality * cp_vap + (1.0 - steam_quality) * cp_liq;

            cp
        },
        FwdEqnRegion::Region5 => cp_tp_5(t, p),
    }
}


/// returns the isochoric (const vol) heat capacity given temperature and pressure
pub fn cv_ph_eqm(p: Pressure, h: AvailableEnergy) -> SpecificHeatCapacity {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => cv_tp_1(t, p),
        FwdEqnRegion::Region2 => cv_tp_2(t, p),
        FwdEqnRegion::Region3 => cv_tp_3(t, p),
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let cv_liq = cv_tp_1(t_sat, p);
            let cv_vap = cv_tp_2(t_sat, p);

            let cv = steam_quality * cv_vap + (1.0 - steam_quality) * cv_liq;

            cv
        },
        FwdEqnRegion::Region5 => cv_tp_5(t, p),
    }
}




/// returns the speed of sound given temperature and pressure
pub fn w_ph_eqm(p: Pressure, h: AvailableEnergy) -> Velocity {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => w_tp_1(t, p),
        FwdEqnRegion::Region2 => w_tp_2(t, p),
        FwdEqnRegion::Region3 => {
            // near supercritical point, we have to be a little bit more careful
            // need to take into account steam quality too

            // so first, test if we are JUST on the saturation line
            // because that is part of the check
            //
            // or rather, just get the volume first
            // this will make it much easier
            let v = v_ph_eqm(p, h);
            let rho = v.recip();
            w_rho_t_3(rho, t)
        },
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let w_liq = w_tp_1(t_sat, p);
            let w_vap = w_tp_2(t_sat, p);

            let w = steam_quality * w_vap + (1.0 - steam_quality) * w_liq;

            w
        },
        FwdEqnRegion::Region5 => w_tp_5(t, p),
    }
}


/// returns the isentropic exponent 
pub fn kappa_ph_eqm(p: Pressure, h: AvailableEnergy) -> Ratio {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => kappa_tp_1(t, p),
        FwdEqnRegion::Region2 => kappa_tp_2(t, p),
        FwdEqnRegion::Region3 => kappa_tp_3(t, p),
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let kappa_liq = kappa_tp_1(t_sat, p);
            let kappa_vap = kappa_tp_2(t_sat, p);

            let kappa = steam_quality * kappa_vap + (1.0 - steam_quality) * kappa_liq;

            kappa
        },
        FwdEqnRegion::Region5 => kappa_tp_5(t, p),
    }
}

/// returns the isobaric cubic expansion coefficient
pub fn alpha_v_ph_eqm(p: Pressure, h: AvailableEnergy) -> TemperatureCoefficient {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => alpha_v_tp_1(t, p),
        FwdEqnRegion::Region2 => alpha_v_tp_2(t, p),
        FwdEqnRegion::Region3 => {
            // near supercritical point, we have to be a little bit more careful
            // need to take into account steam quality too

            // so first, test if we are JUST on the saturation line
            // because that is part of the check
            //
            // or rather, just get the volume first
            // this will make it much easier
            let v = v_ph_eqm(p, h);
            let rho = v.recip();
            alpha_v_rho_t_3(rho, t)
        },
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let alpha_v_liq = alpha_v_tp_1(t_sat, p);
            let alpha_v_vap = alpha_v_tp_2(t_sat, p);

            let alpha_v = steam_quality * alpha_v_vap + (1.0 - steam_quality) * alpha_v_liq;

            alpha_v
        },
        FwdEqnRegion::Region5 => alpha_v_tp_5(t, p),
    }
}


/// returns the isothermal compressibility
pub fn kappa_t_ph_eqm(p: Pressure, h: AvailableEnergy) -> InversePressure {
    let t = t_ph_eqm(p, h);
    let region = ph_flash_region(p, h);

    match region {
        FwdEqnRegion::Region1 => kappa_t_tp_1(t, p),
        FwdEqnRegion::Region2 => kappa_t_tp_2(t, p),
        FwdEqnRegion::Region3 => kappa_t_tp_3(t, p),
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let kappa_t_liq = kappa_t_tp_1(t_sat, p);
            let kappa_t_vap = kappa_t_tp_2(t_sat, p);

            let kappa_t = steam_quality * kappa_t_vap + (1.0 - steam_quality) * kappa_t_liq;

            kappa_t
        },
        FwdEqnRegion::Region5 => kappa_t_tp_5(t, p),
    }
}

/// obtains steam quality (vap fraction) given 
/// pressure and enthalpy 
pub fn x_ph_flash(p: Pressure, h: AvailableEnergy,) -> f64 {

    let region = ph_flash_region(p, h);
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
            // the enthalpy is below critical enthalpy 
            let h_crit = h_3a3b_backwards_ph_boundary(p);

            if h < h_crit {
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
                let h_liq = h_tp_1(t_sat, p);
                let h_vap = h_tp_2(t_sat, p);

                let x: Ratio = (h-h_liq)/(h_vap-h_liq);
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

                let h_liq = h_rho_t_3(v_liq.recip(), t_sat);
                let h_vap = h_rho_t_3(v_vap.recip(), t_sat);

                if h_liq == h_vap {
                    // at supercritical point, just assume vapour
                    // this prevents a divide by zero error
                    return 1.0;
                };


                let x: Ratio = (h-h_liq)/(h_vap-h_liq);
                return x.get::<ratio>();


            }

        },
        // this is a placeholder, 
        // but technically if in region 5, the vapour 
        // quality is 1.0
        FwdEqnRegion::Region5 => 1.0,
    }

}



// allows the user to check which region one is in based on a ph flash
//
// note that ph flash does not work in region 5
pub fn ph_flash_region(p: Pressure, h: AvailableEnergy) -> FwdEqnRegion {

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
    // regimes this is the two phase region
    // this is from h = 1670.9 kJ/kg to about 2563.6 kJ/kg
    // and from 16.529 MPa to 22.064 Mpa (critical point)
    // here is where things are potentially two phase

    let is_region_1 
        = is_ph_point_region_1_and_above_16_529_mpa(p, h);

    if is_region_1 {
        return FwdEqnRegion::Region1;
    };

    let is_region_2 
        = is_ph_point_region_2_and_above_16_529_mpa(p, h);
    if is_region_2 {
        return FwdEqnRegion::Region2;
    };

    // the checks for if it is region 1 or 2 already exclude points 
    // outside h = 1670.9 kJ/kg to about 2563.6 kJ/kg
    //
    // this is above 22.064 MPa

    let is_region_3_above_supercrit 
        = is_ph_point_region_3_and_above_critical_point(p, h);

    if is_region_3_above_supercrit {
        return FwdEqnRegion::Region3;
    };

    // we want to check if this is region 3 and in the range 
    // 16.520 MPa up to crit point 22.064 MPa
    let is_region_3_below_supercrit
        = is_ph_point_region_3_and_from_16_529_mpa_to_crit_temp(p, h);

    if is_region_3_below_supercrit {
        return FwdEqnRegion::Region3;
    };

    // now we shall have to decide if it is region 3 or 4

    let is_region_4 = 
        is_ph_point_region_4_and_above_16_529_mpa(p, h);
    if is_region_4 {
        return FwdEqnRegion::Region4;
    };

    

    // otherwise it's region 3 

    return FwdEqnRegion::Region3;
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

/// viscosity 
pub use crate::dynamic_viscosity::mu_ph_eqm as mu_ph_eqm;
