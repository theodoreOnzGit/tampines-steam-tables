use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::megapascal;
use uom::si::f64::*;

use crate::region_2_vapour::{h_2b2c, p_2b2c};

#[test]
pub fn boundary_2b2c_pressure(){
    // this is the ph point that must be met
    let p_verification = Pressure::new::<megapascal>(0.100_000_000e3);
    let h_verification = AvailableEnergy::new::<kilojoule_per_kilogram>(0.351_600_432_3e4);

    let p_test = p_2b2c(h_verification);

    approx::assert_relative_eq!(
        p_verification.get::<megapascal>(),
        p_test.get::<megapascal>(),
        max_relative=1e-8
        );
}

#[test]
pub fn boundary_2b2c_enthalpy(){
    // this is the ph point that must be met
    let p_verification = Pressure::new::<megapascal>(0.100_000_000e3);
    let h_verification = AvailableEnergy::new::<kilojoule_per_kilogram>(0.351_600_432_3e4);

    let h_test = h_2b2c(p_verification);

    approx::assert_relative_eq!(
        h_verification.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );

}
