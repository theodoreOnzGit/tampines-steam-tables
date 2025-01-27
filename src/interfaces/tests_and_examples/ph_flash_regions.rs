use uom::si::thermodynamic_temperature::kelvin;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::f64::*;

use crate::interfaces::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::interfaces::functional_programming::pt_flash_eqm::FwdEqnRegion;
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;

#[test] 
pub fn test_pts_outside_validity_range(){

    // this is just a reminder to test points 
    // outside validity range of ph flash 
    // eg region 5
    //
    // probably want to assert a failure or something
    //
    // so maybe i want to use an error or something
    todo!();
}


/// below 16.529 MPa 
/// h = 1000 kJ/kg 
///
/// should be region 1
#[test] 
pub fn region_1_test_1(){

    let reference_region = FwdEqnRegion::Region1;
    let p = Pressure::new::<megapascal>(16.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1000.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}


/// below 27.400 MPa 
/// h = 1500 kJ/kg 
///
/// should be region 1
#[test] 
pub fn region_1_test_2(){

    let reference_region = FwdEqnRegion::Region1;
    let p = Pressure::new::<megapascal>(27.400);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1500.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}

/// saturated liquid line 
/// below 16.529 MPa 
/// should be region 1 
#[test] 
pub fn region_1_test_3_sat_liq_line(){

    let reference_region = FwdEqnRegion::Region1;
    let p = Pressure::new::<megapascal>(10.0);
    let t = sat_temp_4(p);
    let h = h_tp_1(t, p);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}


/// saturated liquid line 
/// below 16.529 MPa 
/// should be region 1 
#[test] 
pub fn region_1_test_4_t_623_15k_isotherm(){

    let reference_region = FwdEqnRegion::Region1;
    let p = Pressure::new::<megapascal>(84.0);
    let t = ThermodynamicTemperature::new::<kelvin>(623.15);
    let h = h_tp_1(t, p);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}
