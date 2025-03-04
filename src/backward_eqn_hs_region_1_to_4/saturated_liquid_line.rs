use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

/// this is for eq 2.40 on page 80
const H1_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 27] = [
    [0.0, 14.0, 0.332_171_191_705_237],
    [0.0, 36.0, 0.611_217_706_323_496e-3],
    [1.0, 3.0, -0.882_092_478_906_822e1],
    [1.0, 16.0, -0.455_628_192_543_250],
    [2.0, 0.0, -0.263_483_840_850_452e-4],
    [2.0, 5.0, -0.223_949_661_148_062e2],
    [3.0, 4.0, -0.428_398_660_164_013e1],
    [3.0, 36.0, -0.616_679_338_856_916],
    [4.0, 4.0, -0.146_823_031_104_040e2],
    [4.0, 16.0, 0.284_523_138_727_299e3],
    [4.0, 24.0, -0.113_398_503_195_444e3],
    [5.0, 18.0, 0.115_671_380_760_859e4],
    [5.0, 24.0, 0.395_551_267_359_325e3],
    [7.0, 1.0, -0.154_891_257_229_285e1],
    [8.0, 4.0, 0.194_486_637_751_291e2],
    [12.0, 2.0, -0.357_915_139_457_043e1],
    [12.0, 4.0, -0.335_369_414_148_819e1],
    [14.0, 1.0, -0.664_426_796_332_460],
    [14.0, 22.0, 0.323_321_885_383_934e5],
    [16.0, 10.0, 0.331_766_744_667_084e4],
    [20.0, 12.0, -0.223_501_257_931_087e5],
    [20.0, 28.0, 0.573_953_875_852_936e7],
    [22.0, 8.0, 0.173_226_193_407_919e3],
    [24.0, 3.0, -0.363_968_822_121_321e-1],
    [28.0, 0.0, 0.834_596_332_878_346e-6],
    [32.0, 6.0, 0.503_611_916_682_674e1],
    [32.0, 8.0, 0.655_444_787_064_505e2],
];

/// this is for eq 2.40 on page 80
const H3A_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 19] = [
    [0.0, 1.0, 0.822_673_364_673_336],
    [0.0, 4.0, 0.181_977_213_534_479],
    [0.0, 10.0, -0.112_000_260_313_624e-1],
    [0.0, 16.0, -0.746_778_287_048_033e-3],
    [2.0, 1.0, -0.179_046_263_257_381],
    [3.0, 36.0, 0.424_220_110_836_657e-1],
    [4.0, 3.0, -0.341_355_823_438_768],
    [4.0, 16.0, -0.209_881_740_853_565e1],
    [5.0, 20.0, -0.822_477_343_323_596e1],
    [5.0, 36.0, -0.499_684_082_076_008e1],
    [6.0, 4.0, 0.191_413_958_471_069],
    [7.0, 2.0, 0.581_062_241_093_136e-1],
    [7.0, 28.0, -0.165_505_498_701_029e4],
    [7.0, 32.0, 0.158_870_443_421_201e4],
    [10.0, 14.0, -0.850_623_535_172_818e2],
    [10.0, 32.0, -0.317_714_386_511_207e5],
    [10.0, 36.0, -0.945_890_406_632_871e5],
    [32.0, 0.0, -0.139_273_847_088_690e-5],
    [32.0, 6.0, 0.631_052_532_240_980],
];

/// this function represents the saturated liquid line
/// for hs flashing between region 1 and region 4
pub fn h1_prime_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in H1_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        eta += ni * (sigma - 1.09).powf(ii) * (sigma + 0.366e-4).powf(ji);
    }

    return h_ref * eta;

}



/// this function represents the saturated liquid line
/// for hs flashing between region 3a and region 4
pub fn h3a_prime_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in H3A_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        eta += ni * (sigma - 1.09).powf(ii) * (sigma + 0.366e-4).powf(ji);
    }

    return h_ref * eta;

}
