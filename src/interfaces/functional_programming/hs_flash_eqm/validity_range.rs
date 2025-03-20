use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::{f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};
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
/// this only applies to region 1
/// 
pub fn hs_is_below_isobar_p_100_mpa_in_region1(
    h: AvailableEnergy, s: SpecificHeatCapacity) -> bool {

    // from fig 2.14, I'm going to use the graph to help 
    // h max is 3715.2 kj/kg

    let h_max = AvailableEnergy::new::<kilojoule_per_kilogram>(1670.9);
    let s_max = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.778);

    if h > h_max || s > s_max {
        // in this case we are outside region 1
        return false;
    };
    todo!();
}


/// critical entropy
pub(crate) fn s_crit() -> SpecificHeatCapacity {
    SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        4.412_021_482_234_76
    )

}
