
use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin};
use validity_range::s_crit;


use crate::{backward_eqn_hs_region_1_to_4::{region_1_and_3::hb13_s_boundary_enthalpy, saturated_liquid_line::{h1_prime_s_boundary_enthalpy, h3a_prime_s_boundary_enthalpy}, saturated_vapour_line::h2c3b_prime_s_boundary_enthalpy}, interfaces::functional_programming::ps_flash_eqm::h_ps_eqm};

use super::pt_flash_eqm::FwdEqnRegion;

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
        hs_region_low_entropy_region_1_and_4(h, s);
    };

    let s_bound_low_entropy_region_1_3a_and_4 = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.778);

    if s <= s_bound_low_entropy_region_1_3a_and_4 {
        hs_region_low_entropy_region_1_3a_and_4(h, s);
    };

    let s_crit = s_crit();

    if s < s_crit {
        hs_region_near_crit_entropy_region_3a_and_4(h, s);
    };

    // boundary line belongs to region 3a
    if s == s_crit {
        hs_region_crit_entropy_region_3a_and_4(h, s);
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
        hs_region_near_crit_entropy_region_3b_and_4(h, s);
    };

    if s <= s_b23_max {
        hs_region_near_crit_entropy_region_2c_3b_and_4(h, s);
    };
    
    let s_2bc = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            5.85
        );

    if s <= s_2bc {
        hs_region_high_entropy_region_2c_and_4(h, s);
    };

    // based on fig 2.19
    let s_where_2b_starts = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            6.070
        );

    if s < s_where_2b_starts {
        hs_region_high_entropy_region_2b_and_4(h, s);
    };

    // based on fig 2.19
    let s_where_2a_starts = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            7.852
        );

    if s <= s_where_2a_starts {
        hs_region_high_entropy_region_2b_2a_and_4(h, s);
    };

    // otherwise, should be in region 2a 
    
    hs_region_high_entropy_region_2a_only(h, s)

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

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677);
    let upper_bound_enthalpy = h_ps_eqm(upper_bound_pressure, s);
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

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677);
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

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677);
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

    if h < region_3a_and_4_h_boundary {
        return BackwdEqnSubRegion::Region3a;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_crit_entropy_region_3a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {


    let upper_bound_pressure = Pressure::new::<megapascal>(100.0);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677);
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

    if h < region_3a_and_4_h_boundary {
        return BackwdEqnSubRegion::Region3a;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_near_crit_entropy_region_3b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    let upper_bound_pressure = Pressure::new::<megapascal>(100.0);
    let lower_bound_pressure = Pressure::new::<megapascal>(0.000_611_212_677);
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
    let h2c3b_boundary = h2c3b_prime_s_boundary_enthalpy(s) - 
        AvailableEnergy::new::<kilojoule_per_kilogram>(0.0073);


    if h < h2c3b_boundary {
        return BackwdEqnSubRegion::Region3b;
    } else {
        return BackwdEqnSubRegion::Region4;
    };

}
fn hs_region_near_crit_entropy_region_2c_3b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    todo!();

}
fn hs_region_high_entropy_region_2c_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    todo!();

}
fn hs_region_high_entropy_region_2b_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    todo!();

}
fn hs_region_high_entropy_region_2b_2a_and_4(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    todo!();

}
fn hs_region_high_entropy_region_2a_only(
    h: AvailableEnergy, s: SpecificHeatCapacity,) -> BackwdEqnSubRegion {

    todo!();

}

/// note:
/// (h,s) flashes along the isotherms 273.15K are not implemented 
/// for simplicity to avoid iterations
pub mod validity_range;
