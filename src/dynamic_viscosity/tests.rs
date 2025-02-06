use uom::si::{f64::*, mass_density::kilogram_per_cubic_meter, thermodynamic_temperature::kelvin};

use crate::dynamic_viscosity::psi_1_viscosity;

use super::psi_0_viscosity;
#[test]
pub fn psi_0_viscosity_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let psi_0_1 = 0.970_904_552e1;
    let psi_0_2 = 0.326_046_811e2;
    let psi_0_3 = 0.244_558_002e2;
    let psi_0_4 = 0.441_936_611e2;

    let psi_test = psi_0_viscosity(t1);
    approx::assert_relative_eq!(
        psi_0_1,
        psi_test,
        max_relative=1e-7
        );
    let psi_test = psi_0_viscosity(t2);
    approx::assert_relative_eq!(
        psi_0_2,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_0_viscosity(t3);
    approx::assert_relative_eq!(
        psi_0_3,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_0_viscosity(t4);
    approx::assert_relative_eq!(
        psi_0_4,
        psi_test,
        max_relative=1e-8
        );

}


#[test]
pub fn psi_1_viscosity_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.549_921_814e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.612_391_201e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_847e2);

    let psi_1_1 = 0.916_694_207e2;
    let psi_1_2 = 0.104_200_938e1;
    let psi_1_3 = 0.296_900_349e1;
    let psi_1_4 = 0.102_422_324e1;

    let psi_test = psi_1_viscosity(t1, rho1);
    approx::assert_relative_eq!(
        psi_1_1,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_1_viscosity(t2,rho2);
    approx::assert_relative_eq!(
        psi_1_2,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_1_viscosity(t3,rho3);
    approx::assert_relative_eq!(
        psi_1_3,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_1_viscosity(t4,rho4);
    approx::assert_relative_eq!(
        psi_1_4,
        psi_test,
        max_relative=1e-8
        );

}
