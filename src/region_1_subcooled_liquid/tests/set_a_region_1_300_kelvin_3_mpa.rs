use uom::si::pressure::megapascal;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;

use crate::region_1_subcooled_liquid::v_tp_1;

#[test] 
pub fn specific_vol_regression_set_a(){
    let ref_vol_m3_per_kg = 0.100215168e-2;
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let p = Pressure::new::<megapascal>(3.0);

    let specific_vol_test_m3_per_kg = 
        v_tp_1(t, p).get::<cubic_meter_per_kilogram>();

    approx::assert_relative_eq!(
        ref_vol_m3_per_kg,
        specific_vol_test_m3_per_kg,
        max_relative=1e-9);

    
}
