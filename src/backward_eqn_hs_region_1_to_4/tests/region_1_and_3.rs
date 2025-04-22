use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use crate::backward_eqn_hs_region_1_to_4::region_1_and_3::hb13_s_boundary_enthalpy;

#[test] 
pub fn hs_boundary_b13_prime_eq_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1.632_525_047e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.7);

    let h_test = hb13_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_b13_prime_eq_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1.593_027_214e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.6);

    let h_test = hb13_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_b13_prime_eq_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1.566_104_611e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.5);

    let h_test = hb13_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}



