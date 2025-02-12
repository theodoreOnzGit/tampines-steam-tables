use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::thermodynamic_temperature::kelvin;

use crate::region_1_subcooled_liquid::t_ps_1;


#[test] 
pub fn t_ph_flash_3_mpa_0_5_kj_per_kg_k_entropy(){
    let p = Pressure::new::<megapascal>(3.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.5);
    let t_ref_kelvin = 0.307_842_258e3;

    let t_calc_kelvin = t_ps_1(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_ref_kelvin,
        t_calc_kelvin,
        max_relative=1e-8
        );
}

#[test] 
pub fn t_ph_flash_80_mpa_0_5_kj_per_kg_k_entropy(){
    let p = Pressure::new::<megapascal>(80.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.5);
    let t_ref_kelvin = 0.309_979_785e3;

    let t_calc_kelvin = t_ps_1(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_ref_kelvin,
        t_calc_kelvin,
        max_relative=1e-8
        );
}

#[test] 
pub fn t_ph_flash_80_mpa_3_kj_per_kg_k_entropy(){
    let p = Pressure::new::<megapascal>(80.0);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.0);
    let t_ref_kelvin = 0.565_899_909e3;

    let t_calc_kelvin = t_ps_1(p, s)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_ref_kelvin,
        t_calc_kelvin,
        max_relative=1e-8
        );
}

