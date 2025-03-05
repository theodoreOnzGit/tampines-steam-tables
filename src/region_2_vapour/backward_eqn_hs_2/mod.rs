use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::pressure::megapascal;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

#[inline]
pub fn p_hs_2(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure {

    // assuming is region 2, 
    // if s is less than s2bc, then it's region 2c 
    // this is on page 91 of book
    // the points on this line belong to region 2b
    let s2b2c = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.85);

    if s < s2b2c {
        return p_hs_2c(h, s);
    };

    // next, we use the boundary for h, this is the isobar P = 4 MPa 
    // if a point falls on this boundary line, then it belongs to 
    // subregion 2a
    let h2a2b = h_2a2b(s);

    if h > h2a2b {
        return p_hs_2b(h, s);
    } else {
        return p_hs_2a(h, s);
    };

}

/// boundary equations for 2a and 2b
pub mod boundary_eqns;
pub use boundary_eqns::*;

#[inline]
pub(crate) fn p_hs_2a(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure{

    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(4200.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(12.0);
    let p_ref = Pressure::new::<megapascal>(4.0);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut pi: f64 = 0.0;

    for coeffs in SUBREGION_2A_BACK_COEFFS_HS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (eta - 0.5).powf(ii) * (sigma - 1.2).powf(ji);
    };

    return pi.powi(4) * p_ref;

}

#[inline]
pub(crate) fn p_hs_2b(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure{

    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(4100.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.9);
    let p_ref = Pressure::new::<megapascal>(100.0);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut pi: f64 = 0.0;

    for coeffs in SUBREGION_2B_BACK_COEFFS_HS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (eta - 0.6).powf(ii) * (sigma - 1.01).powf(ji);
    };

    return pi.powi(4) * p_ref;

}

#[inline]
pub(crate) fn p_hs_2c(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure{

    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(3500.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.9);
    let p_ref = Pressure::new::<megapascal>(100.0);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut pi: f64 = 0.0;

    for coeffs in SUBREGION_2C_BACK_COEFFS_HS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (eta - 0.7).powf(ii) * (sigma - 1.1).powf(ji);
    };

    return pi.powi(4) * p_ref;

}

/// based on table 2.82
const SUBREGION_2A_BACK_COEFFS_HS: [[f64; 3]; 29] = [
    [0.0, 1.0, -0.182_575_361_923_032e-1],
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


/// based on table 2.83
const SUBREGION_2B_BACK_COEFFS_HS: [[f64; 3]; 33] = [
    [0.0, 0.0, 0.801_496_989_929_459e-1],
    [0.0, 1.0, -0.543_862_807_146_111],
    [0.0, 2.0, 0.337_455_597_421_283],
    [0.0, 4.0, 0.890_555_451_157_450e1],
    [0.0, 8.0, 0.313_840_736_431_485e3],
    [1.0, 0.0, 0.797_367_065_977_789],
    [1.0, 1.0, -0.121_616_973_556_240e1],
    [1.0, 2.0, 0.872_803_386_937_477e1],
    [1.0, 3.0, -0.169_769_781_757_602e2],
    [1.0, 5.0, -0.186_552_827_328_146e3],
    [1.0, 12.0, 0.951_159_274_344_237e5],
    [2.0, 1.0, -0.189_168_510_120_494e2],
    [2.0, 6.0, -0.433_407_037_194_840e4],
    [2.0, 18.0, 0.543_212_633_012_715e9],
    [3.0, 0.0, 0.144_793_408_386_013],
    [3.0, 1.0, 0.128_024_559637_516e3],
    [3.0, 7.0, -0.672_309_534_071_268e5],
    [3.0, 12.0, 0.336_972_380_095_287e8],
    [4.0, 1.0, -0.586_634_196_762_720e3],
    [4.0, 16.0, -0.221_403_224_769_889e11],
    [5.0, 1.0, 0.171_606_668_708_389e4],
    [5.0, 12.0, -0.570_817_595_806_302e9],
    [6.0, 1.0, -0.312_109_693_178_482e4],
    [6.0, 8.0, -0.207_841_384_633_010e7],
    [6.0, 18.0, 0.305_605_946_157_786e13],
    [7.0, 1.0, 0.322_157_004_314_333e4],
    [7.0, 16.0, 0.326_810_259_797_295e12],
    [8.0, 1.0, -0.144_104_158_934_487e4],
    [8.0, 3.0, 0.410_694_867_802_691e3],
    [8.0, 14.0, 0.109_077_066_873_024e12],
    [8.0, 18.0, -0.247_964_654_258_893e14],
    [12.0, 10.0, 0.188_801_906_865_134e10],
    [14.0, 16.0, -0.123_651_009_018_773e15],
];


/// based on table 2.84
const SUBREGION_2C_BACK_COEFFS_HS: [[f64; 3]; 31] = [
    [0.0, 0.0, 0.112_225_607_199_012],
    [0.0, 1.0, -0.339_005_953_606_712e1],
    [0.0, 2.0, -0.320_503_911_730_094e2],
    [0.0, 3.0, -0.197_597_305_104_900e3],
    [0.0, 4.0, -0.407_693_861_553_446e3],
    [0.0, 8.0, 0.132_943_775_222_331e5],
    [1.0, 0.0, 0.170_846_839_774_007e1],
    [1.0, 2.0, 0.373_694_198_142_245e2],
    [1.0, 5.0, 0.358_144_365_815_434e4],
    [1.0, 8.0, 0.423_014_446_424_664e6],
    [1.0, 14.0, -0.751_071_025_760_063e9],
    [2.0, 2.0, 0.523_446_127_607_898e2],
    [2.0, 3.0, -0.228_351_290_812_417e3],
    [2.0, 7.0, -0.960_652_417_056_937e6],
    [2.0, 10.0, -0.807_059_292_526_074e8],
    [2.0, 18.0, 0.162_698_017_225_669e13],
    [3.0, 0.0, 0.772_465_073_604_171],
    [3.0, 5.0, 0.463_929_973_837_746e5],
    [3.0, 8.0, -0.137_317_885_134_128e8],
    [3.0, 16.0, 0.170_470_392_630_512e13],
    [3.0, 18.0, -0.251_104_628_187_308e14],
    [4.0, 18.0, 0.317_748_830_835_520e14],
    [5.0, 1.0, 0.538_685_623_675_312e2],
    [5.0, 4.0, -0.553_089_094_625_169e5],
    [5.0, 6.0, -0.102_861_522_421_405e7],
    [5.0, 14.0, 0.204_249_418_756_234e13],
    [6.0, 8.0, 0.273_918_446_626_997e9],
    [6.0, 18.0, -0.263_963_146_312_685e16],
    [10.0, 7.0, -0.107_890_854_108_088e10],
    [12.0, 7.0, -0.296_492_620_980_124e11],
    [16.0, 10.0, -0.111_754_907_323_424e16],
];
