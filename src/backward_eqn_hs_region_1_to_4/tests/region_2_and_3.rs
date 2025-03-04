use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use crate::backward_eqn_hs_region_1_to_4::region_2_and_3::tb23_prime_s_boundary_enthalpy;

#[test] 
pub fn hs_boundary_b23_prime_eq_1(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.1);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(7.135_259_364e2);

    let t_test = tb23_prime_s_boundary_enthalpy(s, h);

    approx::assert_relative_eq!(
        t_ref.get::<kelvin>(),
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_b23_prime_eq_2(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.15);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(7.685_345_532e2);

    let t_test = tb23_prime_s_boundary_enthalpy(s, h);

    approx::assert_relative_eq!(
        t_ref.get::<kelvin>(),
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );
}

#[test] 
pub fn hs_boundary_b23_prime_eq_3(){
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.2);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(8.176_202_120e2);

    let t_test = tb23_prime_s_boundary_enthalpy(s, h);

    approx::assert_relative_eq!(
        t_ref.get::<kelvin>(),
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );
}




