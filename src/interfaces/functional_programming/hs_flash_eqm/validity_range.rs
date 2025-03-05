use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::temperature_interval::millikelvin;

use crate::region_1_subcooled_liquid::{p_hs_1, t_ph_1};

/// based on page 72 boundary, we use this
/// for all pressure
pub fn hs_is_above_isotherm_t_273_15_kelvin(
    h: AvailableEnergy, s: SpecificHeatCapacity) -> bool {

    // todo: cover the range of specific enthalpy and specific entropy

    let t_ref = ThermodynamicTemperature::new::<kelvin>(273.15);

    // basically take t1(p1(h,s),h) and check against 273.15
    let p_test_region_1 = p_hs_1(h, s);

    let t1_test = t_ph_1(p_test_region_1, h);

    // note that we add 24 millikelvin to this
    let delta_t_correction = TemperatureInterval::new::<millikelvin>(24.0);

    if t1_test + delta_t_correction >= t_ref {
        return true;
    } else {
        return false;
    };

}


/// based on page 73 boundary, we use this at 
/// t = 273.15 K to t = 1073.15 K
pub fn hs_is_above_isobar_p_100_mpa(
    h: AvailableEnergy, s: SpecificHeatCapacity) -> bool {
    todo!();
}
