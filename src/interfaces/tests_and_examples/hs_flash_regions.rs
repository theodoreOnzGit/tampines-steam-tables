use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};

use crate::interfaces::functional_programming::hs_flash_eqm::{hs_flash_region, BackwdEqnSubRegion};


#[test]
fn region1_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.02);

    let correct_region = BackwdEqnSubRegion::Region1;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


#[test]
fn region1_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1550.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.398);

    let correct_region = BackwdEqnSubRegion::Region1;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}




#[test]
fn region4_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1000.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.398);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


#[test]
fn region4_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1553.9);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.778);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near critical entropy
#[test]
fn region4_test_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1553.9);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.412);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near region 3b 2c boundary
#[test]
fn region4_test_4(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.211);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 2b boundary
#[test]
fn region4_test_5(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.040);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2b and 2a boundary
#[test]
fn region4_test_6(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.0);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near pure 2a region
#[test]
fn region4_test_7(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2500.9);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(8.0);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}

// near critical point, 3a
#[test]
fn region3a_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2087.5);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.4);

    let correct_region = BackwdEqnSubRegion::Region3a;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 3a and 1 boundary
#[test]
fn region3a_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.9);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.778);

    let correct_region = BackwdEqnSubRegion::Region3a;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near crit pt
#[test]
fn region3b_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2300.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.5);

    let correct_region = BackwdEqnSubRegion::Region3b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
#[test]
fn region3b_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2563.6);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.210);

    let correct_region = BackwdEqnSubRegion::Region3b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
#[test]
fn region3b_test_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2563.6);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.048);

    let correct_region = BackwdEqnSubRegion::Region3b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
// grey box top left
#[test]
fn region3b_test_4(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2750.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.07);

    let correct_region = BackwdEqnSubRegion::Region3b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
// grey box bottom left
#[test]
fn region3b_test_5(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.10);

    let correct_region = BackwdEqnSubRegion::Region3b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
// grey box top right
#[test]
fn region2c_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.10);

    let correct_region = BackwdEqnSubRegion::Region2c;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
// grey box top right
#[test]
fn region2c_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.20);

    let correct_region = BackwdEqnSubRegion::Region2c;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
// grey box bottom right
#[test]
fn region2c_test_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.26);

    let correct_region = BackwdEqnSubRegion::Region2c;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
// grey box top left bound
#[test]
fn region2c_test_4(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2812.9);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.261);

    let correct_region = BackwdEqnSubRegion::Region2c;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 3b boundary 
#[test]
fn region2c_test_5(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2812.9);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.3);

    let correct_region = BackwdEqnSubRegion::Region2c;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near 2c 2b boundary 
#[test]
fn region2b_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.86);

    let correct_region = BackwdEqnSubRegion::Region2b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near isotherm 1073.15 K and isobar 100 MPa
#[test]
fn region2b_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3715.2);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.040);

    let correct_region = BackwdEqnSubRegion::Region2b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// near isotherm 1073.15 K and isobar 100 MPa
#[test]
fn region2b_test_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3715.2);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(6.1);

    let correct_region = BackwdEqnSubRegion::Region2b;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// just a random region in 2a
#[test]
fn region2a_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3715.2);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(8.1);

    let correct_region = BackwdEqnSubRegion::Region2a;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}


// just a at this entropy, 
// the x=1 saturation line meets the p=ps line
// but it's in 2a
#[test]
fn region2a_test_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3715.2);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.156);

    let correct_region = BackwdEqnSubRegion::Region2a;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}
