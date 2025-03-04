use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use crate::backward_eqn_hs_region_1_to_4::saturated_vapour_line::h2c3b_prime_s_boundary_enthalpy;
use crate::backward_eqn_hs_region_1_to_4::saturated_vapour_line::h2ab_prime_s_boundary_enthalpy;

#[test] 
pub fn hs_boundary_2a2b_prime_eq_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2.723_729_985e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.0);

    let h_test = h2ab_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_2a2b_prime_eq_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2.599_047_210e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(8.0);

    let h_test = h2ab_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_2a2b_prime_eq_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2.511_861_477e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.0);

    let h_test = h2ab_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}


#[test] 
pub fn hs_boundary_2c3b_prime_eq_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2.687_693_850e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.5);

    let h_test = h2c3b_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_2c3b_prime_eq_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2.451_623_609e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.0);

    let h_test = h2c3b_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_2c3b_prime_eq_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2.144_360_448e3);

    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.5);

    let h_test = h2c3b_prime_s_boundary_enthalpy(s);

    approx::assert_relative_eq!(
        h.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );
}

