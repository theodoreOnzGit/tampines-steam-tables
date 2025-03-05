use uom::si::thermodynamic_temperature::kelvin;
use uom::si::pressure::{megapascal, pascal};
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

#[inline]
pub fn t_hs_2(h: AvailableEnergy, s: SpecificHeatCapacity) -> ThermodynamicTemperature {

    todo!()
}

/// boundary equations for 2a and 2b
pub mod boundary_eqns;
pub use boundary_eqns::*;
