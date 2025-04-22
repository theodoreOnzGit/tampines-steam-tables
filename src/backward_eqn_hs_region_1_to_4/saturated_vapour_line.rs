use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

/// this is for eq 2.40 on page 80
const H2AB_DOUBLE_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 30] = [
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
    [0.0, 0.0, 0.104_351_280_732_769e1],
    [0.0, 3.0, -0.227_807_912_708_513e1],
    [0.0, 4.0, 0.180_535_256_723_202e1],
    [1.0, 0.0, 0.420_440_834_792_042],
    [1.0, 12.0, -0.105_721_244_834_660e6],
    [5.0, 36.0, 0.436_911_607_493_884e25],
    [6.0, 12.0, -0.328_032_702_839_753e12],
    [7.0, 16.0, -0.678_686_760_804_270e16],
    [8.0, 2.0, 0.743_957_464_645_363e4],
    [8.0, 20.0, -0.356_896_445_355_761e20],
    [12.0, 32.0, 0.167_590_585_186_801e32],
    [16.0, 36.0, -0.355_028_625_419_105e38],
    [22.0, 2.0, 0.396_611_982_166_538e12],
    [22.0, 32.0, -0.414_716_268_484_468e41],
    [24.0, 7.0, 0.359_080_103_867_382e19],
    [36.0, 20.0, -0.116_994_334_851_995e41],
];

/// this function represents the saturated liquid line
/// for hs flashing between region 2a and 2b
pub fn h2ab_double_prime_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref_1 = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.21);
    let s_ref_2 = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.2);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let sigma_1: f64 = (s/s_ref_1).get::<ratio>();
    let sigma_2: f64 = (s/s_ref_2).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in H2AB_DOUBLE_PRIME_S_BOUNDARY_EQN_COEFFS {
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

