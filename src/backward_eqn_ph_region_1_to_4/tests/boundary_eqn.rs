use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal};

use crate::backward_eqn_ph_region_1_to_4::p_s3_h;
#[test]
pub fn boundary_eqn_ps3_test_1(){
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let p_ref = Pressure::new::<megapascal>(1.724_175_718_e1);

    let p_test = p_s3_h(h_ref);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn boundary_eqn_ps3_test_2(){
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2000.0);
    let p_ref = Pressure::new::<megapascal>(2.193442957e1);

    let p_test = p_s3_h(h_ref);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
#[test]
pub fn boundary_eqn_ps3_test_3(){
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2400.0);
    let p_ref = Pressure::new::<megapascal>(2.018090839e1);

    let p_test = p_s3_h(h_ref);

    approx::assert_relative_eq!(
        p_ref.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );

}
