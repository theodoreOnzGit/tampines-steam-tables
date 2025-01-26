use uom::si::f64::*;
use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::thermodynamic_temperature::kelvin;

use crate::constants;


/// Returns the region-3 delta
/// dimensionless specific density
/// Specific mass is assumed to be in kg/m3
pub fn delta_3(rho: MassDensity) -> f64 {
    rho.get::<kilogram_per_cubic_meter>() / constants::RHO_C_KG_PER_M3
}
/// Returns the region-3 tau
/// dimensionless temperature
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
pub fn tau_3(t: ThermodynamicTemperature) -> f64 {
    constants::T_C_KELVIN / t.get::<kelvin>()
}

