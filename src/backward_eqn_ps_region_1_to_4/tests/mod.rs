use uom::si::{f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin};

use crate::backward_eqn_ps_region_1_to_4::boundary_eqn_ps3::p_s3_s;

#[test]
pub fn boundary_eqn_ps3_test_1(){
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let p_ref = Pressure::new::<megapascal>(1.687_755_057e1);

    let p_test = p_s3_s(s_ref);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn boundary_eqn_ps3_test_2(){
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.2);
    let p_ref = Pressure::new::<megapascal>(2.164_451_789e1);

    let p_test = p_s3_s(s_ref);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn boundary_eqn_ps3_test_3(){
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.2);
    let p_ref = Pressure::new::<megapascal>(1.668_968_482e1);

    let p_test = p_s3_s(s_ref);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
