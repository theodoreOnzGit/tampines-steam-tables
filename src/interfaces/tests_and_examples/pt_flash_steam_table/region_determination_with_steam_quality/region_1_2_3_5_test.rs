use crate::{constants::p_crit_water, region_4_vap_liq_equilibrium::sat_pressure_4};

pub use super::*;

#[test]
pub fn region1_test(){
    let reference_region = FwdEqnRegion::Region1;

    let t = ThermodynamicTemperature::new::<kelvin>(473.15);

    // below crit point, region 1 
    let p_sat = sat_pressure_4(t);
    let p1 = p_sat + Pressure::new::<bar>(5.0);
    let p_crit = p_crit_water();
    let p2 = p_crit + Pressure::new::<bar>(50.0);

    let steam_quality = 0.0;
    // all these should return region 1
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p1, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p_crit, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p2, steam_quality)
        );
}

#[test]
pub fn region2_test(){
    let reference_region = FwdEqnRegion::Region2;

    let t = ThermodynamicTemperature::new::<kelvin>(873.15);

    // below crit point, region 1 
    let p_sat = sat_pressure_4(t);
    let p1 = p_sat + Pressure::new::<bar>(5.0);
    let p_crit = p_crit_water();
    let p2 = p_crit + Pressure::new::<bar>(50.0);

    let steam_quality = 0.0;
    // all these should return region 1
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p1, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p_crit, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p2, steam_quality)
        );
}


#[test]
pub fn region3_test(){
    let reference_region = FwdEqnRegion::Region3;

    let t = ThermodynamicTemperature::new::<kelvin>(633.15);

    // below crit point, region 1 
    let p_sat = sat_pressure_4(t);
    let p1 = p_sat + Pressure::new::<bar>(0.05);
    let p_crit = p_crit_water();
    let p2 = p_crit + Pressure::new::<bar>(50.0);

    let steam_quality = 0.0;
    // all these should return region 1
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p1, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p_crit, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p2, steam_quality)
        );
}


#[test]
pub fn region5_test(){
    let reference_region = FwdEqnRegion::Region5;

    let t = ThermodynamicTemperature::new::<kelvin>(1193.15);

    // below crit point, region 1 
    let p_sat = sat_pressure_4(t);
    let p1 = p_sat + Pressure::new::<bar>(0.05);
    let p_crit = p_crit_water();
    let p2 = p_crit + Pressure::new::<bar>(10.0);

    let steam_quality = 0.0;
    // all these should return region 1
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p1, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p_crit, steam_quality)
        );
    assert_eq!(
        reference_region,
        region_fwd_eqn_two_phase(t, p2, steam_quality)
        );
}
