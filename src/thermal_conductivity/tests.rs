use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;

use super::*;

#[test]
pub fn lambda_0_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let lambda_0_1 = 0.184_341_883e2;
    let lambda_0_2 = 0.791_034_659e2;
    let lambda_0_3 = 0.545_433_367e2;
    let lambda_0_4 = 0.119_586_108e3;

    let lambda_0_test = lambda_0(t1);
    approx::assert_relative_eq!(
        lambda_0_1,
        lambda_0_test,
        max_relative=1e-8
        );
    let lambda_0_test = lambda_0(t2);
    approx::assert_relative_eq!(
        lambda_0_2,
        lambda_0_test,
        max_relative=1e-8
        );
    let lambda_0_test = lambda_0(t3);
    approx::assert_relative_eq!(
        lambda_0_3,
        lambda_0_test,
        max_relative=1e-8
        );
    let lambda_0_test = lambda_0(t4);
    approx::assert_relative_eq!(
        lambda_0_4,
        lambda_0_test,
        max_relative=1e-8
        );
    
}
#[test]
pub fn lambda_1_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.260_569_558e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.523_371_289e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_848e2);

    let lambda_1_1 = 0.329_016_833e2;
    let lambda_1_2 = 0.110_043_337e1;
    let lambda_1_3 = 0.726_398_725e1;
    let lambda_1_4 = 0.115_280_540e1;

    let lambda_1_test = lambda_1(t1,rho1);
    approx::assert_relative_eq!(
        lambda_1_1,
        lambda_1_test,
        max_relative=1e-8
        );
    let lambda_1_test = lambda_1(t2,rho2);
    approx::assert_relative_eq!(
        lambda_1_2,
        lambda_1_test,
        max_relative=1e-8
        );
    let lambda_1_test = lambda_1(t3,rho3);
    approx::assert_relative_eq!(
        lambda_1_3,
        lambda_1_test,
        max_relative=1e-8
        );
    let lambda_1_test = lambda_1(t4,rho4);
    approx::assert_relative_eq!(
        lambda_1_4,
        lambda_1_test,
        max_relative=1e-8
        );
    
}
