
use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin, thermodynamic_temperature::kelvin};
use validity_range::s_crit;


use crate::{backward_eqn_hs_region_1_to_4::{region_1_and_3::hb13_s_boundary_enthalpy, region_2_and_3::tb23_s_boundary_enthalpy, saturated_liquid_line::{h1_prime_s_boundary_enthalpy, h3a_prime_s_boundary_enthalpy}, saturated_vapour_line::{h2ab_double_prime_s_boundary_enthalpy, h2c3b_prime_s_boundary_enthalpy}}, interfaces::functional_programming::ps_flash_eqm::h_ps_eqm, region_2_vapour::{h_2a2b, p_hs_2c}, region_3_single_phase_plus_supercritical_steam::p_boundary_2_3, region_4_vap_liq_equilibrium::sat_pressure_4};

use super::pt_flash_eqm::{h_tp_eqm_single_phase, FwdEqnRegion};

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
    let upper_bound_enthalpy = h_tp_eqm_single_phase(t_max, p_sat);
    let lower_bound_enthalpy = h_ps_eqm(lower_bound_pressure, s);

    if h > upper_bound_enthalpy {
        panic!("enthalpy too high for hs flash");
    };

    if h < lower_bound_enthalpy {
        dbg!(&(h,lower_bound_enthalpy));
        panic!("enthalpy too low for hs flash");
    };

    let h2ab_double_prime_boundary = h2ab_double_prime_s_boundary_enthalpy(s);

    if h <= h2ab_double_prime_boundary {
        return BackwdEqnSubRegion::Region4;
    } else {
        return BackwdEqnSubRegion::Region2a;
    };



}

/// note:
/// (h,s) flashes along the isotherms 273.15K are not implemented 
/// for simplicity to avoid iterations
pub mod validity_range;

#[cfg(test)]
mod region_split_tests;
