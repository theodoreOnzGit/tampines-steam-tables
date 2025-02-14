use uom::si::{f64::*, force::newton, length::meter, radiant_exposure::joule_per_square_meter, thermodynamic_temperature::kelvin};

use super::water_surf_tension;

/// tests according to table 3.9
#[test]
fn surf_tension_unit_test_1(){
    let t = ThermodynamicTemperature::new::<kelvin>(300.0);
    let sigma_ref = Force::new::<newton>(0.716_859_625e-1)/
        Length::new::<meter>(1.0);

    let sigma_test = water_surf_tension(t);


    approx::assert_relative_eq!(
        sigma_ref.get::<joule_per_square_meter>(),
        sigma_test.get::<joule_per_square_meter>(),
        max_relative=1e-8
        );


}

/// tests according to table 3.9
#[test]
fn surf_tension_unit_test_2(){
    let t = ThermodynamicTemperature::new::<kelvin>(450.0);
    let sigma_ref = Force::new::<newton>(0.428_914_992e-1)/
        Length::new::<meter>(1.0);

    let sigma_test = water_surf_tension(t);


    approx::assert_relative_eq!(
        sigma_ref.get::<joule_per_square_meter>(),
        sigma_test.get::<joule_per_square_meter>(),
        max_relative=1e-8
        );


}


/// tests according to table 3.9
#[test]
fn surf_tension_unit_test_3(){
    let t = ThermodynamicTemperature::new::<kelvin>(600.0);
    let sigma_ref = Force::new::<newton>(0.837_561_087e-2)/
        Length::new::<meter>(1.0);

    let sigma_test = water_surf_tension(t);


    approx::assert_relative_eq!(
        sigma_ref.get::<joule_per_square_meter>(),
        sigma_test.get::<joule_per_square_meter>(),
        max_relative=1e-8
        );


}
