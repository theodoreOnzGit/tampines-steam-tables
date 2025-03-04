use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

/// this is for eq 2.40 on page 80
const H2AB_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 30] = [
    [1.0, 8.0, -0.524_581_170_928_788e3],
    [1.0, 24.0, -0.926_947_218_142_218e7],
    [2.0, 4.0, -0.237_385_107_491_666e3],
    [2.0, 32.0, 0.210_770_155_812_776e11],
    [4.0, 1.0, -0.239_494_562_010_986e2],
    [4.0, 2.0, 0.221_802_480_294_197e3],
    [7.0, 7.0, -0.510_472_533_393_438e7],
    [8.0, 5.0, 0.124_981_396_109_147e7],
    [8.0, 12.0, 0.200_008_436_996_201e10],
    [10.0, 1.0, -0.815_158_509_791_035e3],
    [12.0, 0.0, -0.157_612_685_637_523e3],
    [12.0, 7.0, -0.114_200_422_332_791e11],
    [18.0, 10.0, 0.662_364_680_776_872e16],
    [20.0, 12.0, -0.227_622_818_296_144e19],
    [24.0, 32.0, -0.171_048_081_348_406e32],
    [28.0, 8.0, 0.660_788_766_938_091e16],
    [28.0, 12.0, 0.166_320_055_886_021e23],
    [28.0, 20.0, -0.218_003_784_381_501e30],
    [28.0, 22.0, -0.787_276_140_295_618e30],
    [28.0, 24.0, 0.151_062_329_700_346e32],
    [32.0, 2.0, 0.795_732_170_300_541e7],
    [32.0, 7.0, 0.131_957_647_355_347e16],
    [32.0, 12.0, -0.325_097_068_299_140e24],
    [32.0, 14.0, -0.418_600_611_419_248e26],
    [32.0, 24.0, 0.297_478_906_557_467e35],
    [36.0, 10.0, -0.953_588_761_745_473e20],
    [36.0, 12.0, 0.166_957_699_620_939e25],
    [36.0, 20.0, -0.175_407_764_869_978e33],
    [36.0, 22.0, 0.347_581_490_626_396e35],
    [36.0, 28.0, -0.710_971_318_427_851e39],
];

/// this is for eq 2.40 on page 80
const H2C3B_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 16] = [
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
];

/// this function represents the saturated liquid line
/// for hs flashing between region 2a and 2b
pub fn h2ab_prime_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref_1 = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.21);
    let s_ref_2 = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.2);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let sigma_1: f64 = (s/s_ref_1).get::<ratio>();
    let sigma_2: f64 = (s/s_ref_2).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in H2AB_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        eta += ni * (sigma_1.recip() - 0.513).powf(ii) * (sigma_2 - 0.524).powf(ji);
    }

    return h_ref * eta.exp();

}



/// this function represents the saturated liquid line
/// for hs flashing between region 2c and 3b
pub fn h2c3b_prime_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.9);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in H2C3B_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        eta += ni * (sigma - 1.02).powf(ii) * (sigma - 0.726).powf(ji);
    }

    return h_ref * eta.powi(4);

}

