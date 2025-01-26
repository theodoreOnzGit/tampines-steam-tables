use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;

use crate::region_2_vapour::t_ph_2;

pub mod boundary_tests_2b2c;


#[test]
pub fn subregion_2a_t_ph_test_1(){
    let p = Pressure::new::<megapascal>(0.001);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);
    let t_expected_kelvin = 0.534433241e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2a_t_ph_test_2(){
    let p = Pressure::new::<megapascal>(3.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);
    let t_expected_kelvin = 0.575373370e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}
#[test]
pub fn subregion_2a_t_ph_test_3(){
    let p = Pressure::new::<megapascal>(3.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(4000.0);
    let t_expected_kelvin = 0.101077577e4;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}

#[test]
pub fn subregion_2b_t_ph_test_1(){
    let p = Pressure::new::<megapascal>(5.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3500.0);
    let t_expected_kelvin = 0.801299102e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}


#[test]
pub fn subregion_2b_t_ph_test_2(){
    let p = Pressure::new::<megapascal>(5.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(4000.0);
    let t_expected_kelvin = 0.101531583e4;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}

#[test]
pub fn subregion_2b_t_ph_test_3(){
    let p = Pressure::new::<megapascal>(25.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3500.0);
    let t_expected_kelvin = 0.875279054e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );
}


#[test]
pub fn subregion_2c_t_ph_test_1(){
    let p = Pressure::new::<megapascal>(40.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let t_expected_kelvin = 0.743056411e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );
}


#[test]
pub fn subregion_2c_t_ph_test_2(){
    let p = Pressure::new::<megapascal>(60.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(2700.0);
    let t_expected_kelvin = 0.791137067e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );
}
#[test]
pub fn subregion_2c_t_ph_test_3(){
    let p = Pressure::new::<megapascal>(60.0);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(3200.0);
    let t_expected_kelvin = 0.882756860e3;

    let t_test = t_ph_2(p, h);

    approx::assert_relative_eq!(
        t_expected_kelvin,
        t_test.get::<kelvin>(),
        max_relative=1e-8
        );

}
