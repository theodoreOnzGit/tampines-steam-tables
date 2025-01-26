use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal};

use crate::region_3_vap_liq_mixture_plus_critical_point::h_3a3b_backwards_ph_boundary;

#[test]
pub fn ph_3a_3b_boundary(){
    let p_ref = Pressure::new::<megapascal>(25.0);
    let h_ref_kj_per_kg = 2.095_936_454e3;

    let h_test_kj_per_kg 
        = h_3a3b_backwards_ph_boundary(p_ref)
        .get::<kilojoule_per_kilogram>();

    approx::assert_relative_eq!(
        h_ref_kj_per_kg,
        h_test_kj_per_kg,
        max_relative=1e-8
        );
}
