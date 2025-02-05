use uom::si::{pressure::megapascal, ratio::ratio};

use crate::region_3_single_phase_plus_supercritical_steam::{alpha_v_rho_t_3, cp_rho_t_3, cv_rho_t_3, h_rho_t_3, kappa_rho_t_3, kappa_t_rho_t_3, s_rho_t_3, u_rho_t_3, v_tp_3c, v_tp_3r, v_tp_3s, v_tp_3t, v_tp_3u, v_tp_3x, v_tp_3y, v_tp_3z, w_rho_t_3};
use crate::constants::{p_crit_water, t_crit_water, T_C_KELVIN};

use super::*;

pub fn region_fwd_eqn_two_phase(
    t: ThermodynamicTemperature,
    p: Pressure,
    steam_quality: f64) -> FwdEqnRegion {

    // for the two phase version, I need to determine if it's in the 
    // saturation region 
    // otherwise, revert back to the normal boundary 
    //
    // if in the saturation region, then x 
    // should be between 0 and 1

    let mut x = steam_quality;

    if x < 0.0 {
        x = 0.0;
    };
    if x > 1.0 {
        x = 1.0;
    };

    // now that that's settled, let's check the 
    // points
    //
    // if beyond critical point, don't even need to consider two phase
    //
    // this will cover region 5 automatically
    if t >= t_crit_water() {
        return region_fwd_eqn_single_phase(t, p);
    };

    if p >= p_crit_water() {
        return region_fwd_eqn_single_phase(t, p);
    };

    // if below critical temperature and pressure
    // check if it lies within the sat temperature and pressure 

    let p_sat_reg4 = sat_pressure_4(t);

    // if pressure is at exactly the saturation pressure,
    // and steam quality is between 0 and 1
    // then return region 4
    
    let multiphase_steam: bool = x > 0.0 && x < 1.0;
    // for this, I need to give some leeway due to numerical precision error
    //let steam_pressure_equal_sat_pressure = p == p_sat_reg4;
    // at least to within 0.01%
    let steam_pressure_equal_sat_pressure 
        = ((p/p_sat_reg4).get::<ratio>() - 1.0).abs() < 5e-4;



    let temp_above_623_15_k = t.get::<kelvin>() > 623.15;

    if steam_pressure_equal_sat_pressure && multiphase_steam {
        return FwdEqnRegion::Region4;
    };
    // next, we also need to be extra careful about bubble and dew point 
    let dew_point: bool = x == 1.0;
    let bubble_point: bool = x == 0.0;

    // first case, bubble/dew point at or below 623.15
    if !temp_above_623_15_k && steam_pressure_equal_sat_pressure 
        && bubble_point {
            return FwdEqnRegion::Region1;
    };
    if !temp_above_623_15_k && steam_pressure_equal_sat_pressure 
        && dew_point {
            return FwdEqnRegion::Region2;
    };

    // if bubble or dew point is in region 3,
    if temp_above_623_15_k && steam_pressure_equal_sat_pressure 
        && bubble_point {
            return FwdEqnRegion::Region3;
    };
    if temp_above_623_15_k && steam_pressure_equal_sat_pressure 
        && dew_point {
            return FwdEqnRegion::Region3;
    };
    
    // if none of the conditions are satisfied, return single phase
    
    return region_fwd_eqn_single_phase(t, p);

}


/// returns the enthalpy given temperature and pressure
/// single phase only!
pub fn h_tp_eqm_two_phase(
    t: ThermodynamicTemperature, p: Pressure,
    x: f64) -> AvailableEnergy {
    let region = region_fwd_eqn_two_phase(t, p, x);



    // note that if x = 1.0  or 0.0 exactly, 
    // if goes straight to region 3
    // in that case we must be extra careful


    match region {
        FwdEqnRegion::Region1 => h_tp_1(t, p),
        FwdEqnRegion::Region2 => h_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else if t_kelvin <= 646.483 {
                    v_tp_3x(t, p)
                } else {
                    v_tp_3z(t, p)
                }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if t_kelvin <= 646.599 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let h_liq = h_rho_t_3(v_liq.recip(), t);
            let h_vap = h_rho_t_3(v_vap.recip(), t);
            let h_saturated = steam_quality * h_vap + (1.0 - steam_quality) * h_liq;

            let near_critical_point: bool = (t_kelvin - 647.096).abs() < 0.025;
            let near_saturation_line: bool = (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 5e-4 ;
            // if we are ON the saturated line, we must be mindful
            //dbg!(&(h_liq,h_vap,h_saturated,steam_quality));

            dbg!(&near_critical_point);
            dbg!(&(near_saturation_line,h_liq,h_vap,h_saturated));
            if near_critical_point && near_saturation_line {
                // this was intended to get the enthalpy of vapourisation 
                // to zero near critical point
                return h_tp_3(t, p);
            } else if near_saturation_line {
                return h_saturated;
            };


            // otherwise
            h_tp_3(t, p)
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
               let h_liq = h_tp_1(t, p);
                let h_vap = h_tp_2(t, p);

                let h = steam_quality * h_vap + (1.0 - steam_quality) * h_liq;

                h
            } else if t_sat_kelvin < T_C_KELVIN {
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let h_liq = h_rho_t_3(v_liq.recip(), t);
                let h_vap = h_rho_t_3(v_vap.recip(), t);
                
                let h = steam_quality * h_vap + (1.0 - steam_quality) * h_liq;
                h
            } else {
                // this is exactly at crit point
                h_tp_3(t, p)
            }

        },
        FwdEqnRegion::Region5 => h_tp_5(t, p),
    }
}

/// returns the internal energy given temperature and pressure
pub fn u_tp_eqm_two_phase(
    t: ThermodynamicTemperature, p: Pressure,
    x: f64) -> AvailableEnergy {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => u_tp_1(t, p),
        FwdEqnRegion::Region2 => u_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else if t_kelvin <= 646.483 {
                    v_tp_3x(t, p)
                } else {
                    v_tp_3z(t, p)
                }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if t_kelvin <= 646.599 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let u_liq = u_rho_t_3(v_liq.recip(), t);
            let u_vap = u_rho_t_3(v_vap.recip(), t);
            let u_saturated = steam_quality * u_vap + (1.0 - steam_quality) * u_liq;

            let near_critical_point: bool = ((u_vap-u_liq)/u_liq ).get::<ratio>() < 9e-3;
            let near_saturation_line: bool = (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 ;
            // if we are ON the saturated line, we must be mindful
            if near_critical_point && near_saturation_line {
                return u_tp_3(t, p);
            } else if near_saturation_line {
                return u_saturated;
            };


            // otherwise
            u_tp_3(t, p)

        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let u_liq = u_tp_1(t, p);
                let u_vap = u_tp_2(t, p);

                let u = steam_quality * u_vap + (1.0 - steam_quality) * u_liq;

                u
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

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
pub fn s_tp_eqm_two_phase(t: ThermodynamicTemperature, p: Pressure,
    x: f64) -> SpecificHeatCapacity {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => s_tp_1(t, p),
        FwdEqnRegion::Region2 => s_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let s_liq = s_rho_t_3(v_liq.recip(), t);
            let s_vap = s_rho_t_3(v_vap.recip(), t);
            let s_saturated = steam_quality * s_vap + (1.0 - steam_quality) * s_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 3e-4 {
                return s_saturated;
            };

            // otherwise
            s_tp_3(t, p)
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let s_liq = s_tp_1(t, p);
                let s_vap = s_tp_2(t, p);

                let s = steam_quality * s_vap + (1.0 - steam_quality) * s_liq;

                s
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
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
pub fn cp_tp_eqm_two_phase(t: ThermodynamicTemperature, 
    p: Pressure, x: f64) -> SpecificHeatCapacity {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => cp_tp_1(t, p),
        FwdEqnRegion::Region2 => cp_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let cp_liq = cp_rho_t_3(v_liq.recip(), t);
            let cp_vap = cp_rho_t_3(v_vap.recip(), t);
            let cp_saturated = steam_quality * cp_vap + (1.0 - steam_quality) * cp_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 {
                return cp_saturated;
            };

            // otherwise
            cp_tp_3(t, p)
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let cp_liq = cp_tp_1(t, p);
                let cp_vap = cp_tp_2(t, p);

                let cp = steam_quality * cp_vap + (1.0 - steam_quality) * cp_liq;

                cp
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let cp_liq = cp_rho_t_3(v_liq.recip(), t);
                let cp_vap = cp_rho_t_3(v_vap.recip(), t);
                let cp = steam_quality * cp_vap + (1.0 - steam_quality) * cp_liq;

                cp
            }

        },
        FwdEqnRegion::Region5 => cp_tp_5(t, p),
    }
}


/// returns the isochoric (const vol) heat capacity given temperature and pressure
pub fn cv_tp_eqm_two_phase(t: ThermodynamicTemperature, p: Pressure,
    x: f64) -> SpecificHeatCapacity {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => cv_tp_1(t, p),
        FwdEqnRegion::Region2 => cv_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let cv_liq = cv_rho_t_3(v_liq.recip(), t);
            let cv_vap = cv_rho_t_3(v_vap.recip(), t);
            let cv_saturated = steam_quality * cv_vap + (1.0 - steam_quality) * cv_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 {
                return cv_saturated;
            };

            // otherwise
            cv_tp_3(t, p)
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let cv_liq = cv_tp_1(t, p);
                let cv_vap = cv_tp_2(t, p);

                let cv = steam_quality * cv_vap + (1.0 - steam_quality) * cv_liq;

                cv
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let cv_liq = cv_rho_t_3(v_liq.recip(), t);
                let cv_vap = cv_rho_t_3(v_vap.recip(), t);
                let cv = steam_quality * cv_vap + (1.0 - steam_quality) * cv_liq;

                cv
            }

        },
        FwdEqnRegion::Region5 => cv_tp_5(t, p),
    }
}




/// returns the specific volume given temperature and pressure
pub fn v_tp_eqm_two_phase(t: ThermodynamicTemperature, 
    p: Pressure,
    x: f64) -> SpecificVolume {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => v_tp_1(t, p),
        FwdEqnRegion::Region2 => v_tp_2(t, p),
        // note that for region 3, the backward eqn is used
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else if t_kelvin <= 646.483 {
                    v_tp_3x(t, p)
                } else {
                    v_tp_3z(t, p)
                }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if t_kelvin <= 646.599 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let v_saturated = steam_quality * v_vap + (1.0 - steam_quality) * v_liq;
            let near_saturation_line = 
                (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 5e-4;

            // extremely near critical point, about 373.707 degc
            // or anything more than 370K degc
            // ph flashing algorithm tends to do better 
            let v_region_3 = v_tp_3(t, p);
            // if not at critical point but 
            // we are ON the saturated line, we must be mindful
            if near_saturation_line {
                return v_saturated;
            };
            // otherwise
            v_region_3
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let v_liq = v_tp_1(t, p);
                let v_vap = v_tp_2(t, p);

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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let v = steam_quality * v_vap + (1.0 - steam_quality) * v_liq;

                v
            }

        },
        FwdEqnRegion::Region5 => v_tp_5(t, p),
    }
}


/// returns the speed of sound given temperature and pressure
///
/// note, for region 4, this is estimated using weighted average of 
/// liquid and vapour phases (ACCURACY NOT GUARANTEED)
pub fn w_tp_eqm_two_phase(t: ThermodynamicTemperature, p: Pressure,
    x: f64) -> Velocity {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => w_tp_1(t, p),
        FwdEqnRegion::Region2 => w_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let w_liq = w_rho_t_3(v_liq.recip(), t);
            let w_vap = w_rho_t_3(v_vap.recip(), t);
            let w_saturated = steam_quality * w_vap + (1.0 - steam_quality) * w_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 {
                return w_saturated;
            };

            // otherwise
            w_tp_3(t, p)
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let w_liq = w_tp_1(t, p);
                let w_vap = w_tp_2(t, p);

                let h = steam_quality * w_vap + (1.0 - steam_quality) * w_liq;

                h
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let w_liq = w_rho_t_3(v_liq.recip(), t);
                let w_vap = w_rho_t_3(v_vap.recip(), t);
                let w = steam_quality * w_vap + (1.0 - steam_quality) * w_liq;

                w
            }

        },
        FwdEqnRegion::Region5 => w_tp_5(t, p),
    }
}


/// returns the isentropic exponent 
pub fn kappa_tp_eqm_two_phase(t: ThermodynamicTemperature, 
    p: Pressure, x: f64) -> Ratio {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => kappa_tp_1(t, p),
        FwdEqnRegion::Region2 => kappa_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let kappa_liq = kappa_rho_t_3(v_liq.recip(), t);
            let kappa_vap = kappa_rho_t_3(v_vap.recip(), t);
            let kappa_saturated = steam_quality * kappa_vap + (1.0 - steam_quality) * kappa_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 {
                return Ratio::new::<ratio>(kappa_saturated);
            };

            // otherwise
            kappa_tp_3(t, p)

        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let kappa_liq = kappa_tp_1(t, p);
                let kappa_vap = kappa_tp_2(t, p);

                let kappa = steam_quality * kappa_vap + (1.0 - steam_quality) * kappa_liq;

                kappa
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let kappa_liq = kappa_rho_t_3(v_liq.recip(), t);
                let kappa_vap = kappa_rho_t_3(v_vap.recip(), t);
                let kappa = steam_quality * kappa_vap + (1.0 - steam_quality) * kappa_liq;

                Ratio::new::<ratio>(kappa)
            }

        },
        FwdEqnRegion::Region5 => kappa_tp_5(t, p),
    }
}

/// returns the isobaric cubic expansion coefficient
pub fn alpha_v_tp_eqm_two_phase(t: ThermodynamicTemperature, 
    p: Pressure, x: f64) -> TemperatureCoefficient {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => alpha_v_tp_1(t, p),
        FwdEqnRegion::Region2 => alpha_v_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let alpha_v_liq = alpha_v_rho_t_3(v_liq.recip(), t);
            let alpha_v_vap = alpha_v_rho_t_3(v_vap.recip(), t);
            let alpha_v_saturated = steam_quality * alpha_v_vap + (1.0 - steam_quality) * alpha_v_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 {
                return alpha_v_saturated;
            };

            // otherwise
            alpha_v_tp_3(t, p)
        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let alpha_v_liq = alpha_v_tp_1(t, p);
                let alpha_v_vap = alpha_v_tp_2(t, p);

                let alpha_v = steam_quality * alpha_v_vap + (1.0 - steam_quality) * alpha_v_liq;

                alpha_v
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let alpha_v_liq = alpha_v_rho_t_3(v_liq.recip(), t);
                let alpha_v_vap = alpha_v_rho_t_3(v_vap.recip(), t);
                let alpha_v = steam_quality * alpha_v_vap + (1.0 - steam_quality) * alpha_v_liq;

                alpha_v
            }

        },
        FwdEqnRegion::Region5 => alpha_v_tp_5(t, p),
    }
}


/// returns the isothermal compressibility
pub fn kappa_t_tp_eqm(t: ThermodynamicTemperature, p: Pressure, x: f64) -> InversePressure {
    let region = region_fwd_eqn_two_phase(t, p, x);

    match region {
        FwdEqnRegion::Region1 => kappa_t_tp_1(t, p),
        FwdEqnRegion::Region2 => kappa_t_tp_2(t, p),
        FwdEqnRegion::Region3 => {
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
            let t_kelvin = t.get::<kelvin>();
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 
            let v_vap: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 640.691 {
                    v_tp_3t(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3r(t, p)
                } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
            };

            let v_liq: SpecificVolume = {
                // this covers up to tsat at 643.15 K
                if t_kelvin <= 634.659 {
                    v_tp_3c(t, p)
                } else if t_kelvin <= 643.15 {
                    v_tp_3s(t, p)
                } else if p_mpa <= 21.9316 {
                    v_tp_3u(t, p)
                } else {
                    v_tp_3y(t, p)
                }


            };

            let kappa_t_liq = kappa_t_rho_t_3(v_liq.recip(), t);
            let kappa_t_vap = kappa_t_rho_t_3(v_vap.recip(), t);
            let kappa_t_saturated = steam_quality * kappa_t_vap + (1.0 - steam_quality) * kappa_t_liq;

            // if we are ON the saturated line, we must be mindful
            if (p_mpa - sat_pressure_4(t).get::<megapascal>()).abs() < 1e-4 {
                return kappa_t_saturated;
            };

            // otherwise
            kappa_t_tp_3(t, p)

        },
        FwdEqnRegion::Region4 => {
            // the only time we get region 4 is in multiphase steam
            // but I'm having this here anyway cause I kiasu
            let mut steam_quality = x;
            if steam_quality < 0.0 {
                steam_quality = 0.0;
            } else if steam_quality > 1.0 {
                steam_quality = 1.0;
            }; 

            let t_sat_kelvin = t.get::<kelvin>();
            if t_sat_kelvin <= 623.15 {
                let kappa_t_liq = kappa_t_tp_1(t, p);
                let kappa_t_vap = kappa_t_tp_2(t, p);

                let kappa_t = steam_quality * kappa_t_vap + (1.0 - steam_quality) * kappa_t_liq;

                kappa_t
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
                        v_tp_3t(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3r(t, p)
                    } else // this covers pressure from 21.0434 Mpa to crit point 
                    if p_mpa <= 21.9010 {
                        v_tp_3x(t, p)
                    } else {
                        v_tp_3z(t, p)
                    }
                };

                let v_liq: SpecificVolume = {
                    // this covers up to tsat at 643.15 K
                    if t_sat_kelvin <= 634.659 {
                        v_tp_3c(t, p)
                    } else if t_sat_kelvin <= 643.15 {
                        v_tp_3s(t, p)
                    } else if p_mpa <= 21.9316 {
                        v_tp_3u(t, p)
                    } else {
                        v_tp_3y(t, p)
                    }


                };

                let kappa_t_liq = kappa_t_rho_t_3(v_liq.recip(), t);
                let kappa_t_vap = kappa_t_rho_t_3(v_vap.recip(), t);
                let kappa_t = steam_quality * kappa_t_vap + (1.0 - steam_quality) * kappa_t_liq;

                kappa_t
            }

        },
        FwdEqnRegion::Region5 => kappa_t_tp_5(t, p),
    }
}

