use uom::si::{f64::*, thermodynamic_temperature::kelvin};

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
