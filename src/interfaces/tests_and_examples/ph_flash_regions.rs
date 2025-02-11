pub

use uom::si::thermodynamic_temperature::kelvin;
use uom::si::{available_energy::kilojoule_per_kilogram, pressure::megapascal};
use uom::si::f64::*;

use crate::interfaces::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::interfaces::functional_programming::pt_flash_eqm::FwdEqnRegion;
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;
use crate::region_3_single_phase_plus_supercritical_steam::t_boundary_2_3;
use crate::region_4_vap_liq_equilibrium::{sat_pressure_4, sat_temp_4};



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


/// saturated liquid line 
/// at 16.529 MPa 
/// should be region 1 
#[test] 
pub fn region_1_test_5_t_623_15k_isotherm_and_sat_liq_line(){

    let reference_region = FwdEqnRegion::Region1;
    let t = ThermodynamicTemperature::new::<kelvin>(623.15);
    let p = sat_pressure_4(t);
    let h = h_tp_1(t, p);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}

// ================== now region 2 tests ===================
///
/// this one is below 16.529 MPa
#[test] 
pub fn region_2_test_1(){

    let reference_region = FwdEqnRegion::Region2;
    let p = Pressure::new::<megapascal>(16.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );

}


///
/// this one is at 16.529 MPa
#[test] 
pub fn region_2_test_2(){

    let reference_region = FwdEqnRegion::Region2;
    let p = Pressure::new::<megapascal>(16.529);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );

}


///
/// this one is above 16.529 MPa
#[test] 
pub fn region_2_test_3(){

    let reference_region = FwdEqnRegion::Region2;
    let p = Pressure::new::<megapascal>(91.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );

}


///
/// this one is above 16.529 MPa
/// at the Tb23 line
#[test] 
pub fn region_2_test_4_tb23_line(){

    let reference_region = FwdEqnRegion::Region2;
    let p = Pressure::new::<megapascal>(91.0);
    let t = t_boundary_2_3(p);
    let h = h_tp_2(t, p);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );

}


///
/// sat vap line below 16.529 MPa
#[test] 
pub fn region_2_test_5_sat_vap_line(){

    let reference_region = FwdEqnRegion::Region2;
    let p = Pressure::new::<megapascal>(16.0);
    let t = sat_temp_4(p);
    let h = h_tp_2(t, p);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );

}


///
/// sat vap line below 16.529 MPa
#[test] 
pub fn region_2_test_6_sat_vap_t_b23(){

    let reference_region = FwdEqnRegion::Region2;
    let t = ThermodynamicTemperature::new::<kelvin>(623.15);
    let p = sat_pressure_4(t);
    let h = h_tp_2(t, p);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );

}


// ================== now region 3 tests ===================
/// region 3 above critical point test
#[test] 
pub fn region_3_test_1(){
    let reference_region = FwdEqnRegion::Region3;
    let p = Pressure::new::<megapascal>(23.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2100.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// region 3 below critical point test
/// but above 16.529 MPa
/// sat liq
///
/// looking at ph graph fig 2.5 roughly
#[test] 
pub fn region_3_test_2(){
    let reference_region = FwdEqnRegion::Region3;
    let p = Pressure::new::<megapascal>(18.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// region 3 below critical point test
/// but above 16.529 MPa
///
/// looking at ph graph fig 2.5 roughly
#[test] 
pub fn region_3_test_3(){
    let reference_region = FwdEqnRegion::Region3;
    let p = Pressure::new::<megapascal>(18.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);

    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// region around h = 1670.9 kJ/kg
/// region around p = 37 MPa
#[test] 
pub fn region_3_test_4(){

    let reference_region = FwdEqnRegion::Region3;
    let p = Pressure::new::<megapascal>(37.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1670.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}


/// region around h = 2563.6 kJ/kg
/// region around p = 37 MPa
#[test] 
pub fn region_3_test_5(){

    let reference_region = FwdEqnRegion::Region3;
    let p = Pressure::new::<megapascal>(37.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );


}

/// region 4
///
/// below 623.15 K VLE
/// below 16.529 Mpa
#[test]
pub fn region_4_test_1(){

    let reference_region = FwdEqnRegion::Region4;
    let p = Pressure::new::<megapascal>(6.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1500.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// region 4
///
/// enthalpy at critical pt
/// below 16.529 Mpa
#[test]
pub fn region_4_test_2(){

    let reference_region = FwdEqnRegion::Region4;
    let p = Pressure::new::<megapascal>(6.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}

/// region 4
///
/// enthalpy just above TB23 line
/// below 16.529 Mpa
#[test]
pub fn region_4_test_3(){

    let reference_region = FwdEqnRegion::Region4;
    let p = Pressure::new::<megapascal>(7.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// region 4
///
/// crit enthalpy
/// pressure just below crit pressure
#[test]
pub fn region_4_test_4(){

    let reference_region = FwdEqnRegion::Region4;
    let p = Pressure::new::<megapascal>(22.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// region 5
///
/// should panic
#[test]
#[should_panic]
pub fn region_5_test_1_should_panic(){

    let reference_region = FwdEqnRegion::Region5;
    let p = Pressure::new::<megapascal>(23.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(4300.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// above 100 mpa
///
/// should panic
#[test]
#[should_panic]
pub fn above_100_mpa_should_panic(){

    let reference_region = FwdEqnRegion::Region5;
    let p = Pressure::new::<megapascal>(121.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(4300.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// enthalpy corresponding to temp below 273.15K 
///
/// should panic
#[test]
#[should_panic]
pub fn above_enthalpy_below_273_15k_should_panic(){

    let reference_region = FwdEqnRegion::Region5;
    let p = Pressure::new::<megapascal>(121.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(-5.0);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


/// pressure below 0.000611 MPa (triple pt pressure) should panic
///
/// should panic
#[test]
#[should_panic]
pub fn above_pressure_below_tripe_pt_pressure_should_panic(){

    let reference_region = FwdEqnRegion::Region5;
    let p = Pressure::new::<megapascal>(6.11e-9);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);


    let test_region = ph_flash_region(p, h);

    assert_eq!(
        reference_region,
        test_region
        );
}


