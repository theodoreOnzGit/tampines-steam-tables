use uom::si::{f64::*, mass_density::kilogram_per_cubic_meter, thermodynamic_temperature::kelvin};

use super::water_dielectric_const_rho_t;
#[test]
fn dielectric_constant_test_1(){
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(0.999_242_866e3);
    let t = ThermodynamicTemperature::new::<kelvin>(298.15);

    let ref_dielectric_const = 0.785_907_250e2;
    let test_dielectric_const = water_dielectric_const_rho_t(rho, t);

    approx::assert_relative_eq!(
        ref_dielectric_const,
        test_dielectric_const,
        max_relative=1e-8
        );

}
#[test]
fn dielectric_constant_test_2(){
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(0.260_569_558e2);
    let t = ThermodynamicTemperature::new::<kelvin>(873.15);

    let ref_dielectric_const = 0.112_620_970e1;
    let test_dielectric_const = water_dielectric_const_rho_t(rho, t);

    approx::assert_relative_eq!(
        ref_dielectric_const,
        test_dielectric_const,
        max_relative=1e-8
        );

}
#[test]
fn dielectric_constant_test_3(){
    let rho = MassDensity::new::<kilogram_per_cubic_meter>(0.523_371_289e3);
    let t = ThermodynamicTemperature::new::<kelvin>(673.15);

    let ref_dielectric_const = 0.103_126_058e2;
    let test_dielectric_const = water_dielectric_const_rho_t(rho, t);

    approx::assert_relative_eq!(
        ref_dielectric_const,
        test_dielectric_const,
        max_relative=1e-8
        );

}
