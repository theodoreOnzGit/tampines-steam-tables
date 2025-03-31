use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};

use super::{hs_flash_region, BackwdEqnSubRegion};

#[test]
fn region2a2b2c_test_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(50.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.02);

    let correct_region = BackwdEqnSubRegion::Region2a;
    let test_region = hs_flash_region(h, s);

    assert_eq!(correct_region, test_region);
    
}
