
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::available_energy::kilojoule_per_kilogram;
use validity_range::s_crit;


use crate::region_1_subcooled_liquid::{p_hs_1, t_ph_1};
use crate::region_4_vap_liq_equilibrium::{sat_pressure_4, tsat_hs_4};
use crate::region_3_single_phase_plus_supercritical_steam::v_ps_flash::v_ps_3b;
use crate::region_3_single_phase_plus_supercritical_steam::v_ps_flash::v_ps_3a;
use crate::region_3_single_phase_plus_supercritical_steam::t_ph_3;
use crate::region_3_single_phase_plus_supercritical_steam::p_hs_3b;
use crate::region_3_single_phase_plus_supercritical_steam::p_hs_3a;
use crate::region_3_single_phase_plus_supercritical_steam::p_boundary_2_3;
use crate::region_2_vapour::{h_2a2b, p_hs_2a, p_hs_2b, p_hs_2c, t_ph_2};
use crate::interfaces::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::x_ph_flash;
use crate::backward_eqn_hs_region_1_to_4::saturated_vapour_line::h2c3b_prime_s_boundary_enthalpy;
use crate::backward_eqn_hs_region_1_to_4::saturated_vapour_line::h2ab_double_prime_s_boundary_enthalpy;
use crate::backward_eqn_hs_region_1_to_4::saturated_liquid_line::h3a_prime_s_boundary_enthalpy;
use crate::backward_eqn_hs_region_1_to_4::saturated_liquid_line::h1_prime_s_boundary_enthalpy;
use crate::backward_eqn_hs_region_1_to_4::region_2_and_3::tb23_s_boundary_enthalpy;
use crate::backward_eqn_hs_region_1_to_4::region_1_and_3::hb13_s_boundary_enthalpy;

use super::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, lambda_ph_eqm, mu_ph_eqm, t_ph_eqm, w_ph_eqm};
use super::pt_flash_eqm::{s_tp_eqm_two_phase, FwdEqnRegion};
use super::pt_flash_eqm::s_tp_eqm_single_phase;
use super::pt_flash_eqm::h_tp_eqm_single_phase;
use super::ps_flash_eqm::v_ps_eqm;

#[derive(Debug,PartialEq, Eq, PartialOrd, Ord)]
/// an enum to help represent the appropriate 
/// regions in the forward equations
pub enum BackwdEqnSubRegion {
    /// this is from T = 273.15 K to T=623.15K 
    /// liquid
    Region1,
    /// this is vapour then line p23/t23 
    /// all the way up to 1073.15 K (800 degC)
    ///
    /// 
    /// higher entropy than 5.85 kJ/(kg K)
    /// but even higher than 2b
    /// includes the boundary line h2ab
    Region2a,
    /// this is vapour then line p23/t23 
    /// all the way up to 1073.15 K (800 degC)
    /// 
    /// higher or equal to entropy than 5.85 kJ/(kg K)
    /// includes boundary line 5.85 kJ/(kg K)
    Region2b,
    /// this is vapour then line p23/t23 
    /// all the way up to 1073.15 K (800 degC)
    ///
    /// lower entropy than 5.85 kJ/(kg K)
    Region2c,
    /// this is supercritical region and 
    /// single phase liquid  / vapour near 
    /// supercritical region
    ///
    /// below or equal to critical entropy
    Region3a,
    /// this is supercritical region and 
    /// single phase liquid  / vapour near 
    /// supercritical region
    /// above critical entropy
    Region3b,
    /// two phase vapour liq equilibrium 
    /// region up to supercritical region
    /// (saturation line, but not including the line itself)
    Region4,
    /// ultra high temperature steam  (more than 800 degC)
    /// 1073.15 K to 2273.15 K 
    /// pressure from triple pt pressure to 500 bar
    Region5,
}

impl Into<FwdEqnRegion> for BackwdEqnSubRegion {
    fn into(self) -> FwdEqnRegion {
        match self {
            BackwdEqnSubRegion::Region1 => FwdEqnRegion::Region1,
            BackwdEqnSubRegion::Region2a => FwdEqnRegion::Region2,
            BackwdEqnSubRegion::Region2b => FwdEqnRegion::Region2,
            BackwdEqnSubRegion::Region2c => FwdEqnRegion::Region2,
            BackwdEqnSubRegion::Region3a => FwdEqnRegion::Region3,
            BackwdEqnSubRegion::Region3b => FwdEqnRegion::Region3,
            BackwdEqnSubRegion::Region4 => FwdEqnRegion::Region4,
            BackwdEqnSubRegion::Region5 => FwdEqnRegion::Region5,
        }
    }
}

/// returns temperature given
/// enthalpy and entropy point
///
///
pub fn t_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> ThermodynamicTemperature {
    let (t,_p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    return t;
}
/// returns pressure given
/// enthalpy and entropy point
///
///
pub fn p_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> Pressure {
    let (_t,p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    return p;
}
/// returns specific volume given
/// enthalpy and entropy point
///
///
pub fn v_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> SpecificVolume {
    let (_t,_p,v,_x) = tpvx_hs_flash_eqm(h, s);

    return v;
}
/// returns quality given
/// enthalpy and entropy point
///
///
pub fn x_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> Ratio {
    let (_t,_p,_v,x) = tpvx_hs_flash_eqm(h, s);

    return x;
}

/// returns cp given 
/// enthalpy and entropy point
/// uses ph flash
pub fn cp_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> SpecificHeatCapacity {
    let (_t,p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    let cp = cp_ph_eqm(p, h);

    return cp;

}
/// returns w (speed of sound) given 
/// enthalpy and entropy point
/// uses ph flash
pub fn w_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> Velocity {
    let (_t,p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    let w = w_ph_eqm(p, h);

    return w;

}

/// returns kappa (isentropic exponent) given 
/// enthalpy and entropy point
/// uses ph flash
pub fn kappa_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> Ratio {
    let (_t,p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    let kappa = kappa_ph_eqm(p, h);

    return kappa;

}
/// returns mu, or sometimes eta (dynamic viscosity) given 
/// enthalpy and entropy point
/// uses ph flash
pub fn mu_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> DynamicViscosity {
    let (_t,p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    let mu = mu_ph_eqm(p, h);

    return mu;

}
/// returns lambda (thermal conductivity) given 
/// enthalpy and entropy point
/// uses ph flash
pub fn lambda_hs_eqm(h: AvailableEnergy, s: SpecificHeatCapacity,) -> ThermalConductivity {
    let (_t,p,_v,_x) = tpvx_hs_flash_eqm(h, s);

    let lambda = lambda_ph_eqm(p, h);

    return lambda;

}


/// returns temperature, pressure, specific volume and quality given 
/// enthalpy and entropy point
///
/// I'm doing this combined function to prevent double calculation
///
#[inline]
pub fn tpvx_hs_flash_eqm(h: AvailableEnergy,
    s: SpecificHeatCapacity,) -> 
(ThermodynamicTemperature, Pressure, SpecificVolume, Ratio) {
    let region = hs_flash_region(h, s);

    match region {
        BackwdEqnSubRegion::Region1 => {
            // page 87 of Kretzchmar textbook
            let pressure =  p_hs_1(h, s); 
            let temperature = t_ph_1(pressure, h);
            // in region 1, we are necessarily liquid,
            // quality is zero
            let quality = Ratio::new::<ratio>(0.0);

            let specific_volume = v_ps_eqm(pressure, s);

            return (temperature, pressure, specific_volume, quality);
        },
        BackwdEqnSubRegion::Region2a => {
            // page 92 to 94 of Kretzchmar textbook
            let pressure =  p_hs_2a(h, s); 
            // 
            let temperature = t_ph_2(pressure, h);
            // in region 2, we are necessarily vapour/gas
            // quality is 1
            let quality = Ratio::new::<ratio>(1.0);
            let specific_volume = v_ps_eqm(pressure, s);

            return (temperature, pressure, specific_volume, quality);
        },
        BackwdEqnSubRegion::Region2b => {
            // page 92 to 94 of Kretzchmar textbook
            let pressure =  p_hs_2b(h, s); 
            // 
            let temperature = t_ph_2(pressure, h);
            // in region 2, we are necessarily vapour/gas
            // quality is 1
            let quality = Ratio::new::<ratio>(1.0);
            let specific_volume = v_ps_eqm(pressure, s);

            return (temperature, pressure, specific_volume, quality);
        },
        BackwdEqnSubRegion::Region2c => {
            // page 92 to 94 of Kretzchmar textbook
            let pressure =  p_hs_2c(h, s); 
            // 
            let temperature = t_ph_2(pressure, h);
            // in region 2, we are necessarily vapour/gas
            // quality is 1
            let quality = Ratio::new::<ratio>(1.0);
            let specific_volume = v_ps_eqm(pressure, s);

            return (temperature, pressure, specific_volume, quality);
        },
        BackwdEqnSubRegion::Region3a => {

            // page 97 onwards of Kretzchmar textbook
            let pressure =  p_hs_3a(h, s); 
            // page 100 onwards of Kretzchmar textbook
            let temperature = t_ph_3(pressure, h);
            // then specific volume  page 99 of Kretzchmar textbook
            let specific_volume = v_ps_3a(pressure, s);

            // quality now 
            let quality = x_ph_flash(pressure, h);
            return (temperature, pressure, specific_volume, quality.into());
        },
        BackwdEqnSubRegion::Region3b => {
            // page 97 onwards of Kretzchmar textbook
            let pressure =  p_hs_3b(h, s); 
            // page 100 onwards of Kretzchmar textbook
            let temperature = t_ph_3(pressure, h);
            // then specific volume  page 99 of Kretzchmar textbook
            let specific_volume = v_ps_3b(pressure, s);
            // quality now 
            let quality = x_ph_flash(pressure, h);
            return (temperature, pressure, specific_volume, quality.into());
        },
        BackwdEqnSubRegion::Region4 => {
            // page 101
            // note, this only works for temperatures below 623.15 K
            // not for temperatures near critical point.
                
            let max_sat_temp_for_backward = 
                ThermodynamicTemperature::new::<kelvin>(623.15);
            let sat_pressure_for_backward = 
                sat_pressure_4(max_sat_temp_for_backward);

            let steam_quality_bound = 1.0;
            let min_entropy_for_backward_eqn = 
                s_tp_eqm_two_phase(max_sat_temp_for_backward, 
                    sat_pressure_for_backward, steam_quality_bound);


            let mut sat_temp = tsat_hs_4(h, s);

            if s >= min_entropy_for_backward_eqn {

                // page 103 
                let sat_pressure = sat_pressure_4(sat_temp);

                // now, we are using the enthalpy, temperature and 
                // pressure to find quality using 
                // x= (h - h_liq)/(h_vap-h_liq)
                // this is the same algorithm I used in the x_ph_flash 
                // calculation
                let quality = x_ph_flash(sat_pressure, h);
                // quality is then used to calculate specific volume...
                // however, it is double calculation here...
                // for x_ph,
                // I'm not overly concerned about computational cost now 
                // but it is an inefficiency
                let specific_volume = v_ps_eqm(sat_pressure, s);
                return (sat_temp, sat_pressure, specific_volume, quality.into());
            } else {

                // if in regime above 623.15 K, 
                // or below the threshold entropy 
                // we need another procedure...
                // may include iteration...
                // to determine temperature
                // page 103 
                let sat_pressure = sat_pressure_4(sat_temp);

                // now, we are using the enthalpy, temperature and 
                // pressure to find quality using 
                // x= (h - h_liq)/(h_vap-h_liq)
                // this is the same algorithm I used in the x_ph_flash 
                // calculation
                let quality = x_ph_flash(sat_pressure, h);
                // quality is then used to calculate specific volume...
                // however, it is double calculation here...
                // for x_ph,
                // I'm not overly concerned about computational cost now 
                // but it is an inefficiency
                let specific_volume = v_ps_eqm(sat_pressure, s);
                return (sat_temp, sat_pressure, specific_volume, quality.into());
            };
        },
        BackwdEqnSubRegion::Region5 => {
            unimplemented!("Region 5 does not have (h,s) flashing");
        },
    }

}

/// allows the user to check which region one is in based on a ph flash
///
/// note that ph flash does not work in region 5
///
/// the way to do region separation is first by entropy according to 
/// fig 2.14
///
/// once that is done, then we separate region by enthalpy.
pub fn hs_flash_region(h: AvailableEnergy, s: SpecificHeatCapacity) -> BackwdEqnSubRegion {

    // this is absolute minimum and max entropy
    // based on fig 2.14
    let s_min = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(-0.00858);
    let s_max = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(11.92);

    if s < s_min {
        panic!("entropy too low for hs flash");
    };
    if s > s_max {
        panic!("entropy too high for hs flash");
    };

    let s_bound_low_entropy_region_1_and_4 = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.398);
    

    if s < s_bound_low_entropy_region_1_and_4 {
        return hs_region_low_entropy_region_1_and_4(h, s);
    };

    let s_bound_low_entropy_region_1_3a_and_4 = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.778);

    if s <= s_bound_low_entropy_region_1_3a_and_4 {
        return hs_region_low_entropy_region_1_3a_and_4(h, s);
    };

    let s_crit = s_crit();

    if s < s_crit {
        return hs_region_near_crit_entropy_region_3a_and_4(h, s);
    };

    // boundary line belongs to region 3a
    if s == s_crit {
        return hs_region_crit_entropy_region_3a_and_4(h, s);
    };

    // this is for the region where TB23 and p2c are used
    let s_b23_min = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.048_096_828
        );

    let s_b23_max = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.260_578_707
        );

    if s < s_b23_min {
        return hs_region_near_crit_entropy_region_3b_and_4(h, s);
    };

    if s <= s_b23_max {
        return hs_region_near_crit_entropy_region_2c_3b_and_4(h, s);
    };
    
    let s_2bc = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.85
        );
    

    if s <= s_2bc {
        return hs_region_high_entropy_region_2c_and_4(h, s);
    };

    // based on fig 2.19
    //
    // also, need to consider 1073.15 K isotherma when s > 6.040 kJ/(kg K)
    let s_where_pure_2b_ends = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            6.070
        );

    if s < s_where_pure_2b_ends {
        return hs_region_high_entropy_region_2b_and_4(h, s);
    };

    // based on fig 2.19
    let s_where_pure_2a_starts = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            7.852
        );

    if s <= s_where_pure_2a_starts {
        return hs_region_high_entropy_region_2b_2a_and_4(h, s);
    };

    // otherwise, should be in region 2a 
    
    return hs_region_high_entropy_region_2a_only(h, s);

}

fn hs_region_low_entropy_region_1_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let region_1_and_4_h_boundary: AvailableEnergy = 
        h1_prime_s_boundary_enthalpy(s)
        - AvailableEnergy::new::<kilojoule_per_kilogram>(0.0034);
    // on page 77 the entire boundary along single phase regions 
    // 1 to 3 and
    // two phase region 4 is considered to belong to both regions 
    //
    // for this calculation though, I'll just assign it to single phase 
    // region by default for simplicity
    //
    // on page 82
    // for points close to boundary, an error of 0.0034 kJ/kg of enthalpy 
    // should be subtracted from the enthalpy to ensure that it is correctly 
    // assigned to single phase region (region 1)

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    dbg!(&upper_bound_enthalpy);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };


    if h >= region_1_and_4_h_boundary {

        return BackwdEqnSubRegion::Region1;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}

fn hs_region_low_entropy_region_1_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };

    // now add corrected h 
    // page 82 
    // helps assign it to single phase region 1
    // rather than two phase region 4

        let region_1_and_4_h_boundary: AvailableEnergy = 
        h1_prime_s_boundary_enthalpy(s)
        - AvailableEnergy::new::<kilojoule_per_kilogram>(0.0034);

    if h < region_1_and_4_h_boundary {
        return BackwdEqnSubRegion::Region4;
    };

    // corrected h for hb13 equation, to assign it to region 1 
    // preferably
    // on page 85

    let hb13_boundary = hb13_s_boundary_enthalpy(s)
        + AvailableEnergy::new::<kilojoule_per_kilogram>(0.018);

    if h > hb13_boundary {
        return BackwdEqnSubRegion::Region3a;
    } else {
        return BackwdEqnSubRegion::Region1;
    };

    

}

fn hs_region_near_crit_entropy_region_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };
    let region_3a_and_4_h_boundary: AvailableEnergy = 
        h3a_prime_s_boundary_enthalpy(s)
        - AvailableEnergy::new::<kilojoule_per_kilogram>(0.0045);

    if h >= region_3a_and_4_h_boundary {
        return BackwdEqnSubRegion::Region3a;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_crit_entropy_region_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {


    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };
    let region_3a_and_4_h_boundary: AvailableEnergy = h3a_prime_s_boundary_enthalpy(s)
        - AvailableEnergy::new::<kilojoule_per_kilogram>(0.0045);

    if h >= region_3a_and_4_h_boundary {
        return BackwdEqnSubRegion::Region3a;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_near_crit_entropy_region_3b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };
    // on page 84 for critical entropy all the way up to 
    // bound of region 2c, we correct with 0.0073 kJ/(kg )
    //
    // but if from scrit to below s''(623.15)K, that is 
    //
    // 4.412 kJ/(kg K) to 5.211 kJ/(kg K) 
    // then correct to 0.0073
    //
    // otherwise  from 5.211 kJ/(kg K) to 5.85 kJ/(kg K)
    // correct to 0.0058
    //
    // for this region 3b and 4, we don't have to worry about discrepancy 
    // in this correction factor yet, just apply 0.0073 kJ/(kg)
    //

    // on page 84 for critical entropy all the way up to 
    // bound of region 2c, we correct with 0.0073 kJ/(kg )
    let h2c3b_boundary = h2c3b_prime_s_boundary_enthalpy(s) - 
        AvailableEnergy::new::<kilojoule_per_kilogram>(0.0073);


    if h >= h2c3b_boundary {
        return BackwdEqnSubRegion::Region3b;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_near_crit_entropy_region_2c_3b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };

    // on page 84 for critical entropy all the way up to 
    // bound of region 2c, we correct with 0.0073 kJ/(kg)
    //
    // but if from scrit to below s''(623.15)K, that is 
    //
    // 4.412 kJ/(kg K) to 5.211 kJ/(kg K) 
    // then correct to 0.0073 kJ/kg
    //
    // otherwise  from 5.211 kJ/(kg K) to 5.85 kJ/(kg K)
    // correct to 0.0058 kJ/kg
    //

    let correction_factor: AvailableEnergy;

    let correction_factor_threshold_entropy = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.211);
    if s < correction_factor_threshold_entropy {
        correction_factor = AvailableEnergy::new::<kilojoule_per_kilogram>(0.0073);
    } else {
        correction_factor = AvailableEnergy::new::<kilojoule_per_kilogram>(0.0058);
    };

    let h2c3b_boundary = 
        h2c3b_prime_s_boundary_enthalpy(s) - 
        correction_factor;

    if h < h2c3b_boundary {
        return BackwdEqnSubRegion::Region4;
    }; 
    // now a few more guard clauses to bound maximum and minimum h 
    let hb23_min = AvailableEnergy::new::<kilojoule_per_kilogram>(
        2.563_592_004e3
    );
    let hb23_max = AvailableEnergy::new::<kilojoule_per_kilogram>(
        2.812_942_061e3
    );
    if h >= hb23_max {
        return BackwdEqnSubRegion::Region2c;
    };
    if h <= hb23_min {
        return BackwdEqnSubRegion::Region3b;
    };


    // now as in page 85 to 86, we first determine tb23(h,s)
    let temp_b23_hs = tb23_s_boundary_enthalpy(s, h);
    // using this temperature, we use pb23 boundary 
    // also, we apply the correction factor on page 80
    let p_b23_boundary = p_boundary_2_3(temp_b23_hs) * (1.0 + 4.5e-5);

    // next we determine the pressure p2c(h,s)

    let p2c_hs = p_hs_2c(h, s);

    if p2c_hs <= p_b23_boundary {
        return BackwdEqnSubRegion::Region2c;
    } else {
        return BackwdEqnSubRegion::Region3b;
    };

}
fn hs_region_high_entropy_region_2c_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };

    // on page 84 for critical entropy all the way up to 
    // bound of region 2c, we correct with 0.0073 kJ/(kg K)
    //
    // but if from scrit to below s''(623.15)K, that is 
    //
    // 4.412 kJ/(kg K) to 5.211 kJ/(kg K) 
    // then correct to 0.0073
    //
    // otherwise  from 5.211 kJ/(kg K) to 5.85 kJ/(kg K)
    // correct to 0.0058
    //
    // here, we don't have to worry about this bound, just correct 
    // using 0.0058 kJ/kg
    let h2c3b_boundary = h2c3b_prime_s_boundary_enthalpy(s) - 
        AvailableEnergy::new::<kilojoule_per_kilogram>(0.0058);

    if h < h2c3b_boundary {
        return BackwdEqnSubRegion::Region4;
    } else {
        return BackwdEqnSubRegion::Region2c;
    };

}
fn hs_region_high_entropy_region_2b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {


    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };

    if s < SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.040) {
        // upper bound enthalpy checking depends on pressure only 
        // when s < 6.040 kJ/(kg K)
        let upper_bound_pressure = Pressure::new::<megapascal>(100.0 - 1.0e-4);
        let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);

        if h > upper_bound_enthalpy {
            panic!("enthalpy too high for hs flash");
        };

    } else {
        // otherwise, the upper bound enthalpy is determined by the 
        // temperature isotherm
        let p_sat = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));
        let t_max = ThermodynamicTemperature::new::<kelvin>(1073.15);
        let upper_bound_enthalpy = h_tp_eqm_single_phase(t_max, p_sat);
        if h > upper_bound_enthalpy {

            dbg!(&(h,upper_bound_enthalpy));
            panic!("enthalpy too high for hs flash");
        };


    };

    let h2ab_double_prime_boundary = h2ab_double_prime_s_boundary_enthalpy(s);

    if h > h2ab_double_prime_boundary {
        return BackwdEqnSubRegion::Region2b;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_high_entropy_region_2b_2a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    // the upper bound enthalpy is determined by the 
    // temperature isotherm
    let p_sat = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));
    let t_max = ThermodynamicTemperature::new::<kelvin>(1073.15);
    let upper_bound_enthalpy = h_tp_eqm_single_phase(t_max, p_sat);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        panic!("enthalpy too low for hs flash");
    };
    // now still need to check if temperature is 
    // too high
    additional_temperature_check_for_1073_15_k_isotherm(h, s);

    let h2ab_double_prime_boundary = h2ab_double_prime_s_boundary_enthalpy(s);

    if h <= h2ab_double_prime_boundary {
        return BackwdEqnSubRegion::Region4;
    }

    let h2ab_boundary = h_2a2b(s);

    // boundary 2ab belongs to subregion 2a 
    // see page 91

    if h <= h2ab_boundary {
        return BackwdEqnSubRegion::Region2a;
    } else {
        return BackwdEqnSubRegion::Region2b;
    };


}
fn hs_region_high_entropy_region_2a_only(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    // the upper bound enthalpy is determined by the 
    // temperature isotherm
    let p_sat = sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));
    let t_max = ThermodynamicTemperature::new::<kelvin>(1073.15);
    let conservative_upper_bound_enthalpy = h_tp_eqm_single_phase(t_max, p_sat);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > conservative_upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        dbg!(&(h,lower_bound_enthalpy));
        panic!("enthalpy too low for hs flash");
    };

    // now still need to check if temperature is 
    // too high
    additional_temperature_check_for_1073_15_k_isotherm(h, s);


    let h2ab_double_prime_boundary = h2ab_double_prime_s_boundary_enthalpy(s);

    if h <= h2ab_double_prime_boundary {
        return BackwdEqnSubRegion::Region4;
    } else {
        return BackwdEqnSubRegion::Region2a;
    };



}


fn additional_temperature_check_for_1073_15_k_isotherm(
    h: AvailableEnergy, s: SpecificHeatCapacity,){

    let t_bound = ThermodynamicTemperature::new::<kelvin>(1073.15);
    let low_bound_entropy_pressure = 
        Pressure::new::<megapascal>(100.0 - 1.0e-3);
    let mid_bound_entropy_pressure = 
        Pressure::new::<megapascal>(4.0);
    let high_bound_entropy_pressure = 
        sat_pressure_4(ThermodynamicTemperature::new::<kelvin>(273.15));
    let low_bound_entropy = 
        s_tp_eqm_single_phase(t_bound, low_bound_entropy_pressure);
    let mid_bound_entropy = 
        s_tp_eqm_single_phase(t_bound, mid_bound_entropy_pressure);
    let high_bound_entropy = 
        s_tp_eqm_single_phase(t_bound, high_bound_entropy_pressure);

    if s < low_bound_entropy {
        panic!("entropy too low for 1073.15K isotherm additional check");
    };
    if s > high_bound_entropy {
        panic!("entropy too high for 1073.15K isotherm additional check");
    };

    if s < mid_bound_entropy {
        let p = p_hs_2b(h, s);
        let t = t_ph_eqm(p, h);
        
        if t > t_bound {
            dbg!(&(t,t_bound));
            dbg!(&(h,s));
            panic!("temperature too high for (h,s) point (2b regime)");
        };
    } else {
        let p = p_hs_2a(h, s);
        let t = t_ph_eqm(p, h);
        
        if t > t_bound {
            dbg!(&(t,t_bound));
            dbg!(&(h,s));
            panic!("temperature too high for (h,s) point (2a regime)");
        };

    };

}

/// note:
/// (h,s) flashes along the isotherms 273.15K are not implemented 
/// for simplicity to avoid iterations
pub mod validity_range;

