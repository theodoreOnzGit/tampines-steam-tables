use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

/// this is for eq 2.44 on page 84
/// based on table 2.73
const HB13_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.913_965_547_600_543],
    [1.0, -2.0, -0.430_944_856_041_991e-4],
    [1.0, 2.0, 0.603_235_694_765_419e2],
    [3.0, -12.0, 0.117_518_273_082_168e-17],
    [5.0, -4.0, 0.220_000_904_781_292],
    [6.0, -3.0, -0.690_815_545_851_641e2],
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
pub fn hb13_prime_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in HB13_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        eta += ni * (sigma - 0.884).powf(ii) * (sigma - 0.864).powf(ji);
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

