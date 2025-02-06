use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;

use super::lambda_0;

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
