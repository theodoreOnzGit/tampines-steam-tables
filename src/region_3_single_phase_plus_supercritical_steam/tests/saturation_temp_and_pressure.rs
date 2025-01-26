use uom::si::pressure::megapascal;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::region_3_single_phase_plus_supercritical_steam::{p_boundary_2_3, t_boundary_2_3};

#[test] 
pub fn sat_temp_to_pressure_p_b23(){

    let ref_pressure_mpa = 0.1652916425e2;
    let t = ThermodynamicTemperature::new::<kelvin>(0.623150000e3);

    let pressure_test_mpa = 
        p_boundary_2_3(t).get::<megapascal>();

    approx::assert_relative_eq!(
        ref_pressure_mpa,
        pressure_test_mpa,
        max_relative=1e-8);
    
}

#[test] 
pub fn sat_pressure_to_temp_t_b23(){

    let pressure_mpa = 0.1652916425e2;
    let ref_t = ThermodynamicTemperature::new::<kelvin>(0.623150000e3);
    let p = Pressure::new::<megapascal>(pressure_mpa);

    let t_test_kelvin = 
        t_boundary_2_3(p).get::<kelvin>();

    approx::assert_relative_eq!(
        ref_t.get::<kelvin>(),
        t_test_kelvin,
        max_relative=1e-8);
    
}

