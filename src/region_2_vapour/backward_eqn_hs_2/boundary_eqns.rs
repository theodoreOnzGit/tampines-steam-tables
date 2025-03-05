use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;


/// eqn for determining pressure boundary between subregion 2b and 2c
/// using dimensionless enthalpy eta
#[inline]
pub fn h_2a2b(s: SpecificHeatCapacity) -> AvailableEnergy {
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        1.0);

    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0);

    let sigma = (s/s_ref).get::<ratio>();


    let n1: f64 = -0.349_898_083_432_139e4;
    let n2: f64 = 0.257_560_716_905_876e4;
    let n3: f64 = -0.421_073_558_227_969e3;
    let n4: f64 = 0.276_349_063_799_944e2;

    let eta: f64 = n1
     + n2 * sigma
     + n3 * sigma.powi(2)
     + n4 * sigma.powi(3);

    return eta * h_ref;
}

#[test]
fn h2a2b_test(){
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(3376.437_884);

    let h_test = h_2a2b(s);

    approx::assert_relative_eq!(
        h_ref.get::<kilojoule_per_kilogram>(),
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-8
        );


}
