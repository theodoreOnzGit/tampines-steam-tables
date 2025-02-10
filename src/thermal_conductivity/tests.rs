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

    let lambda_1_test = lambda_1(rho1,t1);
    approx::assert_relative_eq!(
        lambda_1_1,
        lambda_1_test,
        max_relative=1e-8
        );
    let lambda_1_test = lambda_1(rho2,t2);
    approx::assert_relative_eq!(
        lambda_1_2,
        lambda_1_test,
        max_relative=1e-8
        );
    let lambda_1_test = lambda_1(rho3,t3);
    approx::assert_relative_eq!(
        lambda_1_3,
        lambda_1_test,
        max_relative=1e-8
        );
    let lambda_1_test = lambda_1(rho4,t4);
    approx::assert_relative_eq!(
        lambda_1_4,
        lambda_1_test,
        max_relative=1e-8
        );
    
}

#[test]
fn c_test(){
    let _t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let _t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let _t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let _t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.260_569_558e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.523_371_289e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_848e2);

    let delta_1: f64 = (rho1/rho_crit_water()).get::<ratio>();
    let delta_2: f64 = (rho2/rho_crit_water()).get::<ratio>();
    let delta_3: f64 = (rho3/rho_crit_water()).get::<ratio>();
    let delta_4: f64 = (rho4/rho_crit_water()).get::<ratio>();

    let c_1 = 0.129_592_952e-1;
    let c_2 = 0.163_793_337e0;
    let c_3 = 0.940_881_573e-1;
    let c_4 = 0.168_780_175e0;

    dbg!(&(delta_1,delta_2,delta_3,delta_4));
    let c_test = captial_c(delta_1);
    approx::assert_relative_eq!(
        c_1,
        c_test,
        max_relative=1e-8
        );
    let c_test = captial_c(delta_2);
    approx::assert_relative_eq!(
        c_2,
        c_test,
        max_relative=1e-8
        );
    let c_test = captial_c(delta_3);
    approx::assert_relative_eq!(
        c_3,
        c_test,
        max_relative=1e-8
        );
    let c_test = captial_c(delta_4);
    approx::assert_relative_eq!(
        c_4,
        c_test,
        max_relative=1e-8
        );
    
}


#[test]
fn captial_b_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.260_569_558e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.523_371_289e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_848e2);

    let delta_1: f64 = (rho1/rho_crit_water()).get::<ratio>();
    let delta_2: f64 = (rho2/rho_crit_water()).get::<ratio>();
    let delta_3: f64 = (rho3/rho_crit_water()).get::<ratio>();
    let delta_4: f64 = (rho4/rho_crit_water()).get::<ratio>();

    let theta_1: f64 = (t1/t_crit_water()).get::<ratio>();
    let theta_2: f64 = (t2/t_crit_water()).get::<ratio>();
    let theta_3: f64 = (t3/t_crit_water()).get::<ratio>();
    let theta_4: f64 = (t4/t_crit_water()).get::<ratio>();

    let kappa_t_1 = 0.451_570_597e-3 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_2 = 0.105_138_803e0 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_3 = 0.141_857_631e-1 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_4 = 0.510_625_539e-1 * Pressure::new::<megapascal>(1.0).recip();

    let b_2 = 5.639_822_730e-3;
    let b_3 = 0.373_064_478e0;
    // b1 and b4 values are not given in the table, but they default to 
    // zero because they are outside the range of validity
    let b_1 = 0.0;
    let b_4 = 0.0;

    let n5 = 1.5;

    dbg!(&(delta_1,delta_2,delta_3,delta_4));

    let b_test = captial_b(delta_1, theta_1, kappa_t_1, n5);
    approx::assert_abs_diff_eq!(
        b_1,
        b_test,
        epsilon=0.0
        );
    let b_test = captial_b(delta_2, theta_2, kappa_t_2, n5);
    approx::assert_relative_eq!(
        b_2,
        b_test,
        max_relative=1e-7
        );
    let b_test = captial_b(delta_3, theta_3, kappa_t_3, n5);
    approx::assert_relative_eq!(
        b_3,
        b_test,
        max_relative=1e-8
        );
    let b_test = captial_b(delta_4, theta_4, kappa_t_4, n5);
    approx::assert_abs_diff_eq!(
        b_4,
        b_test,
        epsilon=0.0
        );
    
}


#[test]
fn small_a_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.260_569_558e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.523_371_289e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_848e2);

    let delta_1: f64 = (rho1/rho_crit_water()).get::<ratio>();
    let delta_2: f64 = (rho2/rho_crit_water()).get::<ratio>();
    let delta_3: f64 = (rho3/rho_crit_water()).get::<ratio>();
    let delta_4: f64 = (rho4/rho_crit_water()).get::<ratio>();

    let theta_1: f64 = (t1/t_crit_water()).get::<ratio>();
    let theta_2: f64 = (t2/t_crit_water()).get::<ratio>();
    let theta_3: f64 = (t3/t_crit_water()).get::<ratio>();
    let theta_4: f64 = (t4/t_crit_water()).get::<ratio>();

    let kappa_t_1 = 0.451_570_597e-3 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_2 = 0.105_138_803e0 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_3 = 0.141_857_631e-1 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_4 = 0.510_625_539e-1 * Pressure::new::<megapascal>(1.0).recip();

    // a1 and a4 values are not given in the table, 
    let a_1 = 0.0;
    let a_2 = 0.271_968_296e-1;
    let a_3 = 0.105_363_489e1;
    let a_4 = 0.0;

    let n3 = 0.135_882_142_589_674e1;
    let n4 = 0.508_474_576_271;
    let n5 = 1.5;

    dbg!(&(delta_1,delta_2,delta_3,delta_4));

    let a_test = small_a(n3, delta_1, theta_1, kappa_t_1, n4, n5);
    approx::assert_abs_diff_eq!(
        a_1,
        a_test,
        epsilon=0.0
        );
    let a_test = small_a(n3, delta_2, theta_2, kappa_t_2, n4, n5);
    approx::assert_relative_eq!(
        a_2,
        a_test,
        max_relative=1e-7
        );
    let a_test = small_a(n3, delta_3, theta_3, kappa_t_3, n4, n5);
    approx::assert_relative_eq!(
        a_3,
        a_test,
        max_relative=1e-8
        );
    
    let a_test = small_a(n3, delta_4, theta_4, kappa_t_4, n4, n5);
    approx::assert_abs_diff_eq!(
        a_4,
        a_test,
        epsilon=0.0
        );
}


#[test]
fn capital_a_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.260_569_558e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.523_371_289e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_848e2);

    let delta_1: f64 = (rho1/rho_crit_water()).get::<ratio>();
    let delta_2: f64 = (rho2/rho_crit_water()).get::<ratio>();
    let delta_3: f64 = (rho3/rho_crit_water()).get::<ratio>();
    let delta_4: f64 = (rho4/rho_crit_water()).get::<ratio>();

    let theta_1: f64 = (t1/t_crit_water()).get::<ratio>();
    let theta_2: f64 = (t2/t_crit_water()).get::<ratio>();
    let theta_3: f64 = (t3/t_crit_water()).get::<ratio>();
    let theta_4: f64 = (t4/t_crit_water()).get::<ratio>();

    let kappa_t_1 = 0.451_570_597e-3 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_2 = 0.105_138_803e0 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_3 = 0.141_857_631e-1 * Pressure::new::<megapascal>(1.0).recip();
    let kappa_t_4 = 0.510_625_539e-1 * Pressure::new::<megapascal>(1.0).recip();

    // a1 and a4 values are not given in the table, 
    let a_1 = 0.0;
    let a_2 = 0.917_330_648e-2;
    let a_3 = 0.176_976_803e0;
    let a_4 = 0.0;

    let n2 = 0.636_619_772_367_581;
    let n3 = 0.135_882_142_589_674e1;
    let n4 = 0.508_474_576_271;
    let n5 = 1.5;

    let b1 = 0.101_056_194e1;
    let b2 = 0.133_679_365e1;
    let b3 = 0.294_802_310e1;
    let b4 = 0.128_536_509e1;


    dbg!(&(delta_1,delta_2,delta_3,delta_4));

    let a_test = captial_a(n2, n3, 
        delta_1, theta_1, kappa_t_1, 
        n4, n5, b1);
    approx::assert_abs_diff_eq!(
        a_1,
        a_test,
        epsilon=0.0
        );
    let a_test = captial_a(n2, n3, 
        delta_2, theta_2, kappa_t_2, 
        n4, n5, b2);
    approx::assert_relative_eq!(
        a_2,
        a_test,
        max_relative=1e-7
        );
    let a_test = captial_a(n2, n3, 
        delta_3, theta_3, kappa_t_3, 
        n4, n5, b3);
    approx::assert_relative_eq!(
        a_3,
        a_test,
        max_relative=1e-8
        );
    
    let a_test = captial_a(n2, n3, 
        delta_4, theta_4, kappa_t_4, 
        n4, n5, b4);
    approx::assert_abs_diff_eq!(
        a_4,
        a_test,
        epsilon=0.0
        );
}
