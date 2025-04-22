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
    [-12.0, 2.0, 0.125_244_360_717_979e-12],
    [-12.0, 10.0, -0.126_599_322_553_713e-1],
    [-12.0, 12.0, 0.506_878_030_140_626e1],
    [-12.0, 14.0, 0.317_847_171_154_202e2],
    [-12.0, 20.0, -0.391_041_161_399_932e6],
    [-10.0, 2.0, -0.975_733_406_392_044e-10],
    [-10.0, 10.0, -0.186_312_419_488_279e2],
    [-10.0, 14.0, 0.510_973_543_414_101e3],
    [-10.0, 18.0, 0.373_847_005_822_362e6],
    [-8.0 ,2.0, 0.299_804_024_666_572e-7],
    [-8.0, 8.0, 0.200_544_393_820_342e2],
    [-6.0, 2.0, -0.498_030_487_662_829e-5],
    [-6.0, 6.0, -0.102_301_806_360_030e2],
    [-6.0, 7.0, 0.552_819_126_990_325e2],
    [-6.0, 8.0, -0.206_211_367_510_878e3],
    [-5.0, 10.0, -0.794_012_232_324_823e4],
    [-4.0, 4.0, 0.782_248_472_028_153e1],
    [-4.0, 5.0, -0.586_544_326_902_468e2],
    [-4.0, 8.0, 0.355_073_647_696_481e4],
    [-3.0, 1.0, -0.115_303_107_290_162e-3],
    [-3.0, 3.0, -0.175_092_403_171_802e1],
    [-3.0, 5.0, 0.257_981_687_748_160e3],
    [-3.0, 6.0, -0.727_048_374_179_467e3],
    [-2.0, 0.0, 0.121_644_822_609_198e-3],
    [-2.0, 1.0, 0.393_137_871_762_692e-1],
    [-1.0, 0.0, 0.704_181_005_909_296e-2],
    [0.0, 3.0, -0.829_108_200_698_110e2],
    [2.0, 0.0, -0.265_178_818_131_250],
    [2.0, 1.0, 0.137_531_682_453_991e2],
    [5.0, 0.0, -0.522_394_090_753_046e2],
    [6.0, 1.0, 0.240_556_298_941_048e4],
    [8.0, 1.0, -0.227_361_631_268_929e5],
    [10.0, 1.0, 0.890_746_343_932_567e5],
    [14.0, 3.0, -0.239_234_565_822_486e8],
    [14.0, 7.0, 0.568_795_808_129_714e10],
];
