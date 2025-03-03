use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use crate::region_1_subcooled_liquid::backward_eqn_hs_1::h1_prime_s_boundary_enthalpy;

#[test] 
pub fn hs_boundary_eq_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3.085_509_647e2);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(1.0);

    let h_test = h1_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_eq_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(7.006_304_472e2);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(2.0);

    let h_test = h1_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_eq_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1.198_359_754e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.0);

    let h_test = h1_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}
