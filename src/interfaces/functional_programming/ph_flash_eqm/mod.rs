use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, ratio::ratio, thermodynamic_temperature::kelvin};

use crate::{region_1_subcooled_liquid::{alpha_v_tp_1, cp_tp_1, cv_tp_1, h_tp_1, kappa_t_tp_1, kappa_tp_1, s_tp_1, t_ph_1, u_tp_1, v_tp_1, w_tp_1, InversePressure}, region_2_vapour::{alpha_v_tp_2, cp_tp_2, cv_tp_2, h_tp_2, kappa_t_tp_2, kappa_tp_2, s_tp_2, t_ph_2, u_tp_2, v_tp_2, w_tp_2}, region_3_single_phase_plus_supercritical_steam::{alpha_v_tp_3, cp_tp_3, cv_tp_3, kappa_t_tp_3, kappa_tp_3, s_tp_3, t_ph_3, u_tp_3, v_ph_3, w_tp_3}, region_4_vap_liq_equilibrium::{sat_pressure_4, sat_temp_4}, region_5_superheated_steam::{alpha_v_tp_5, cp_tp_5, cv_tp_5, kappa_t_tp_5, kappa_tp_5, s_tp_5, u_tp_5, w_tp_5}};

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

            let v_liq = v_tp_1(t_sat, p);
            let v_vap = v_tp_2(t_sat, p);

            let v = steam_quality * v_vap + (1.0 - steam_quality) * v_liq;

            v
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
        FwdEqnRegion::Region3 => u_tp_3(t, p),
        FwdEqnRegion::Region4 => {
            // for region 4 specifically, we determine 
            // steam quality first
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let u_liq = u_tp_1(t_sat, p);
            let u_vap = u_tp_2(t_sat, p);

            let u = steam_quality * u_vap + (1.0 - steam_quality) * u_liq;

            u
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
        FwdEqnRegion::Region3 => s_tp_3(t, p),
        FwdEqnRegion::Region4 => {
            // I'm just using quality to interpolate here 
            // not sure if 100% correct
            let steam_quality = x_ph_flash(p, h);
            let t_sat = sat_temp_4(p);

            let s_liq = s_tp_1(t_sat, p);
            let s_vap = s_tp_2(t_sat, p);

            let s = steam_quality * s_vap + (1.0 - steam_quality) * s_liq;

            s
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
        FwdEqnRegion::Region3 => cp_tp_3(t, p),
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
        FwdEqnRegion::Region3 => w_tp_3(t, p),
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
        FwdEqnRegion::Region3 => alpha_v_tp_3(t, p),
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
        FwdEqnRegion::Region2 => 0.0,
        FwdEqnRegion::Region3 => {
            // region 3 is special, if it is equal or above 
            // crit point, then just consider it vapour,
            // doesn't really matter
            if p >= Pressure::new::<megapascal>(22.064) {
                return 1.0;
            };

            // otherwise, in region 3, we check if 
            // the enthalpy is below critical enthalpy 
            let h_crit = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);

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

            let h_liq = h_tp_1(t_sat, p);
            let h_vap = h_tp_2(t_sat, p);

            let x: Ratio = (h-h_liq)/(h_vap-h_liq);

            return x.get::<ratio>();

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
