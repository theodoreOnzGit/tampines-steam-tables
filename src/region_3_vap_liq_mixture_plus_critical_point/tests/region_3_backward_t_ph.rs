use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, specific_volume::cubic_meter_per_kilogram, thermodynamic_temperature::kelvin};

use crate::region_3_vap_liq_mixture_plus_critical_point::{t_ph_3a, t_ph_3b};
#[test] 
pub fn t3a_ph_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let t_verification_kelvin = 6.293083892e2;


    let t_test_kelvin = 
        t_ph_3a(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}

#[test] 
pub fn t3a_ph_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2000.0);
    let t_verification_kelvin = 6.905718338e2;


    let t_test_kelvin = 
        t_ph_3a(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}

#[test] 
pub fn t3a_ph_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2100.0);
    let t_verification_kelvin = 7.336163014e2;


    let t_test_kelvin = 
        t_ph_3a(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}


#[test] 
pub fn t3b_ph_test1(){
    let p = Pressure::new::<megapascal>(20.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2500.0);
    let t_verification_kelvin = 6.418_418_053e2;


    let t_test_kelvin = 
        t_ph_3b(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}

#[test] 
pub fn t3b_ph_test2(){
    let p = Pressure::new::<megapascal>(50.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2400.0);
    let t_verification_kelvin = 7.351_848_618_e2;


    let t_test_kelvin = 
        t_ph_3b(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}

#[test] 
pub fn t3b_ph_test3(){
    let p = Pressure::new::<megapascal>(100.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let t_verification_kelvin = 8.420_460_876_e2;


    let t_test_kelvin = 
        t_ph_3b(p, h)
        .get::<kelvin>();

    approx::assert_relative_eq!(
        t_verification_kelvin,
        t_test_kelvin,
        max_relative=1e-8
        );
}

