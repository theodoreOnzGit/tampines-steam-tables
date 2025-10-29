use uom::si::{f64::*, temperature_interval::kelvin};

// for this dummy DNB type model
//
// I specify a maximum UA, say 700,000 W/K as in the FHR educational simulator 
//
// The minimum UA is perhaps 1,000 W/K. It should start from a minimum value,
// head to maximum at say 
//
// 30 degrees of superheat, then dip back down to 1,000 W/K
//
// We can make a dummy model of bubble accumulation and growth in the wetted 
// portion.
//
// A typical curve could look like this:
//
// https://www.sciencedirect.com/topics/engineering/pool-boiling-curve
// 
// degree of superheat (degc),heat flux (W/m^2),heat transfer coeff (W/m^2 K)
// 1,1.00E+03,1.00E+03
// 5,1.00E+04,2.00E+03
// 12,7.00E+04,5.83E+03
// 30,1.00E+06,3.33E+04
// 320,1.00E+04,3.13E+01
// 2000,1.00E+06,5.00E+02
//
// Planted on log log plots, the heat transfer coeff vs degree of superheat 
// resembles that of a quantum harmonic oscillator (one level above ground 
// state)
// https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
//
// that is x * e^(-a x^2)
//
//
// When plotting log10(degree of superheat) vs 
// log10 (heat transfer coeff)
//
// I make a variable x
//
// x_mod = log10(degree of superheat) - 2 
//
// apparently, this curve fits pretty well
// 
//
// y = -b * x_mod * exp (-a * x_mod * x_mod) + c
//
// These parameters reproduce log10(heat transfer coeff) quite okay. 
// At least fits to within the curve
// a = 2.4
// b = 6 
// c = 3
//
//
// This correlation has a maximum value of around 4.24e4 W/(m^2 K) at 
// 30 degrees kelvin of superheat
//
//
pub(crate) fn pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
    delta_t: TemperatureInterval) -> f64 {

    let degree_of_superheat_kelvin: f64 = 
        delta_t.get::<uom::si::temperature_interval::kelvin>();

    let x_mod: f64 = degree_of_superheat_kelvin.log10() - 2.0;

    let a = 2.4_f64;
    let b = 6.0_f64;
    let c = 3.0_f64;

    let exponent_term = (-a * x_mod * x_mod).exp();

    let y = -b * x_mod * exponent_term + c;

    let heat_transfer_coeff_watt_per_m2_kelvin = 10_f64.powf(y);

    let fraction_of_maximum: f64 = 
        heat_transfer_coeff_watt_per_m2_kelvin/4.24e4;

    return fraction_of_maximum;
}
/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_1(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(1.0);
    let expected_fraction_of_max_htc = 
        0.023629084;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_2(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(5.0);
    let expected_fraction_of_max_htc = 
        0.032133374;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_3(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(12.0);
    let expected_fraction_of_max_htc = 
        0.124350868;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_4(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(30.0);
    let expected_fraction_of_max_htc = 
        1.000846847;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_5(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(100.0);
    let expected_fraction_of_max_htc = 
        0.023584906;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_6(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(320.0);
    let expected_fraction_of_max_htc = 
        0.000536758;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_7(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(1000.0);
    let expected_fraction_of_max_htc = 
        0.006734827;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}


/// this pool boiling function should give roughly these results:
/// degree of superheat(K),fraction of max heat transfer coeff
/// 1,0.023629084
/// 5,0.032133374
/// 12,0.124350868
/// 30,1.000846847
/// 100,0.023584906
/// 320,0.000536758
/// 1000,0.006734827
/// 2000,0.017310593
#[test] 
fn pool_boiling_assert_test_8(){
    let degree_of_superheat_kelvin = 
        TemperatureInterval::new::<kelvin>(2000.0);
    let expected_fraction_of_max_htc = 
        0.017310593;
    let actual_fraction_of_max_htc = 
        pool_boiling_improvised_correlation_as_fraction_of_maximum_lib(
            degree_of_superheat_kelvin
        );

    approx::assert_relative_eq!(
        expected_fraction_of_max_htc,
        actual_fraction_of_max_htc,
        max_relative=1e-7
    );
}
