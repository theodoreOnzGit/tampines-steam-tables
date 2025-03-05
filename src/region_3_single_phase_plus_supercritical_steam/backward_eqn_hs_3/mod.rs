// for region 3 
// backward equations 
// use T(p,h) after finding pressure 
// and v(p,s) after finding pressure
//



use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::available_energy::kilojoule_per_kilogram;

// assuming we are already in region 3
// calculate temperature given p and h
#[inline]
pub fn p_hs_3(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure {

    let s3a3b = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        4.412_021_482_234_76
    );

    // assuming hs point is in region 3, 
    if s > s3a3b {
        return p_hs_3b(h, s);
    } else {
        return p_hs_3a(h, s);
    };

}

pub(crate) fn p_hs_3a(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure {

    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2300.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(4.4);
    let p_ref = Pressure::new::<megapascal>(99.0);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut pi: f64 = 0.0;

    for coeffs in SUBREGION_3A_BACK_COEFFS_HS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (eta - 1.01).powf(ii) * (sigma - 0.75).powf(ji);
    };

    return pi * p_ref;


}


pub(crate) fn p_hs_3b(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure {

    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.3);
    let p_ref = Pressure::new::<megapascal>(16.6);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut pi: f64 = 0.0;

    for coeffs in SUBREGION_3B_BACK_COEFFS_HS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (eta - 0.681).powf(ii) * (sigma - 0.792).powf(ji);
    };

    return pi.recip() * p_ref;



}



/// based on table 2.88
const SUBREGION_3A_BACK_COEFFS_HS: [[f64; 3]; 33] = [
    [0.0, 0.0, 0.770_889_828_326_934e1],
    [0.0, 1.0, -0.260_835_009_128_688e2],
    [0.0, 5.0, 0.267_416_218_930_389e3],
    [1.0, 0.0, 0.172_221_089_496_844e2],
    [1.0, 3.0, -0.293_542_332_145_970e3],
    [1.0, 4.0, 0.614_135_601_882_478e3],
    [1.0, 8.0, -0.610_562_757_725_674e5],
    [1.0, 14.0, -0.651_272_251_118_219e8],
    [2.0, 6.0, 0.735_919_313_521_937e5],
    [2.0, 16.0, -0.116_646_505_914_191e11],
    [3.0, 0.0, 0.355_267_086_434_461e2],
    [3.0, 2.0, -0.596_144_543_825_955e3],
    [3.0, 3.0, -0.475_842_430_145_708e3],
    [4.0, 0.0, 0.696_781_965_359_503e2],
    [4.0, 1.0, 0.335_674_250_377_312e3],
    [4.0, 4.0, 0.250_526_809_130_882e5],
    [4.0, 5.0, 0.146_997_380_630_766e6],
    [5.0, 28.0, 0.538_069_315_091_534e20],
    [6.0, 28.0, 0.143_619_827_291_346e22],
    [7.0, 24.0, 0.364_985_866_165_994e20],
    [8.0, 1.0, -0.254_741_561_156_755e4],
    [10.0, 32.0, 0.240_120_197_096_563e28],
    [10.0, 36.0, -0.393_847_464_679_496e30],
    [14.0, 22.0, 0.147_073_407_024_852e25],
    [18.0, 28.0, -0.426_391_250_432_059e32],
    [20.0, 36.0, 0.194_509_340_621_077e39],
    [22.0, 16.0, 0.666_212_132_114_896e24],
    [22.0, 28.0, 0.706_777_016_552_858e34],
    [24.0, 36.0, 0.175_563_621_975_576e42],
    [28.0, 16.0, 0.108_408_607_429_124e29],
    [28.0, 36.0, 0.730_872_705_175_151e44],
    [32.0, 10.0, 0.159_145_847_398_870e25],
    [32.0, 28.0, 0.377_121_605_943_324e41],
];

/// based on table 2.89
const SUBREGION_3B_BACK_COEFFS_HS: [[f64; 3]; 35] = [
    [0.0, 1.0, -0.182_575_361_923_032e-1],
    [0.0, 1.0, -0.182_575_361_923_032e-1],
    [0.0, 1.0, -0.182_575_361_923_032e-1],
    [0.0, 3.0, -0.125_229_548_799_536],
    [0.0, 3.0, -0.125_229_548_799_536],
    [0.0, 3.0, -0.125_229_548_799_536],
    [0.0, 3.0, -0.125_229_548_799_536],
    [0.0, 3.0, -0.125_229_548_799_536],
    [0.0, 6.0, 0.592_290_437_320_145],
    [0.0, 16.0, 0.604_769_706_185_122e1],
    [0.0, 20.0, 0.238_624_965_444_474e3],
    [0.0, 22.0, -0.298_639_090_222_922e3],
    [1.0, 0.0, 0.512_250_813_040_750e-1],
    [1.0, 1.0, -0.437_266_515_606_486],
    [1.0, 2.0, 0.413_336_902_999_504],
    [1.0, 3.0, -0.516_468_254_574_773e1],
    [1.0, 5.0, -0.557_014_838_445_711e1],
    [1.0, 6.0, 0.128_555_037_824_478e2],
    [1.0, 10.0, 0.114_144_108_953_290e2],
    [1.0, 16.0, -0.119_504_225_652_714e3],
    [1.0, 20.0, -0.284_777_985_961_560e4],
    [1.0, 22.0, 0.431_757_846_408_006e4],
    [2.0, 3.0, 0.112_894_040_802_650e1],
    [2.0, 16.0, 0.197_409_186_206_319e4],
    [2.0, 20.0, 0.151_612_444_706_087e4],
    [3.0, 0.0, 0.141_324_451_421_235e-1],
    [3.0, 2.0, 0.585_501_282_219_601],
    [3.0, 3.0, -0.297_258_075_863_012e1],
    [3.0, 6.0, 0.594_567_314_847_319e1],
    [3.0, 16.0, -0.623_656_565_798_905e4],
    [4.0, 16.0, 0.965_986_235_133_332e4],
    [5.0, 3.0, 0.681_500_934_948_134e1],
    [5.0, 16.0, -0.633_207_286_824_489e4],
    [6.0, 3.0, -0.558_919_224_465_760e1],
    [7.0, 1.0, 0.400_645_798_472_063e-1],
];
