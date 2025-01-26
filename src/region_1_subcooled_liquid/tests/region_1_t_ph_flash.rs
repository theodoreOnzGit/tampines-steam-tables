use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;

use crate::region_1_subcooled_liquid::t_ph_1;

#[test] 
pub fn t_ph_flash_3_mpa_500_kj_per_kg_enthalpy(){
    let p = Pressure::new::<megapascal>(3.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(500.0);
    let t_ref_kelvin = 0.391798509e3;

    let t_calc_kelvin = t_ph_1(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_ref_kelvin,
        t_calc_kelvin,
        max_relative=1e-8
        );
}


#[test] 
pub fn t_ph_flash_80_mpa_500_kj_per_kg_enthalpy(){
    let p = Pressure::new::<megapascal>(80.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(500.0);
    let t_ref_kelvin = 0.378108626e3;

    let t_calc_kelvin = t_ph_1(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_ref_kelvin,
        t_calc_kelvin,
        max_relative=1e-8
        );
}


#[test] 
pub fn t_ph_flash_80_mpa_1500_kj_per_kg_enthalpy(){
    let p = Pressure::new::<megapascal>(80.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1500.0);
    let t_ref_kelvin = 0.611041229e3;

    let t_calc_kelvin = t_ph_1(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_ref_kelvin,
        t_calc_kelvin,
        max_relative=1e-8
        );
}
