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
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2287.5);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(8.0);

    let correct_region = BackwdEqnSubRegion::Region4;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}

