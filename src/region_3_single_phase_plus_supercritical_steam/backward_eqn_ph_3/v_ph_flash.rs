use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, specific_volume::cubic_meter_per_kilogram};

use super::is_3a_when_in_region_3;
/// from table 2.41
const V_PH_SUBREGION_3A_COEFFS: [[f64; 3]; 32] = [
    [-12.0, 6.0, 0.529_944_062_966_028e-2],
    [-12.0, 8.0, -0.170_099_690_234_461],
    [-12.0,12.0, 0.111_323_814_312_927e2],
    [-12.0,18.0, -0.217_898_123_145_125e4],
    [-10.0, 4.0, -0.506_061_827_980_875e-3],
    [-10.0, 7.0, 0.556_495_239_685_324],
    [-10.0,10.0, -0.943_672_726_904_016e1],
    [ -8.0, 5.0, -0.297_856_807_561_527],
    [ -8.0,12.0, 0.939_353_943_717_186e2],
    [ -6.0, 3.0, 0.192_944_939_465_981e-1],
    [ -6.0, 4.0, 0.421_740_664_704_763],
    [ -6.0,22.0, -0.368_914_126_282_330e7],
    [ -4.0, 2.0, -0.737_566_847_600_639e-2],
    [ -4.0, 3.0, -0.354_753_242_424_366],
    [ -3.0, 7.0, -0.199_768_169_338_727e1],
    [ -2.0, 3.0, 0.115_456_297_059_049e1],
    [ -2.0,16.0, 0.568_366_875_815_960e4],
    [ -1.0, 0.0, 0.808_169_540_124_668e-2],
    [ -1.0, 1.0, 0.172_416_341_519_307],
    [ -1.0, 2.0, 0.104_270_175_292_927e1],
    [ -1.0, 3.0, -0.297_691_372_792_847],
    [  0.0, 0.0, 0.560_394_465_163_593],
    [  0.0, 1.0, 0.275_234_661_176_914],
    [  1.0, 0.0, -0.148_347_894_866_012],
    [  1.0, 1.0, -0.651_142_513_478_515e-1],
    [  1.0, 2.0, -0.292_468_715_386_302e1],
    [  2.0, 0.0, 0.664_876_096_952_665e-1],
    [  2.0, 2.0, 0.352_335_014_263_844e1],
    [  3.0, 0.0, -0.146_340_792_313_332e-1],
    [  4.0, 2.0, -0.224_503_486_668_184e1],
    [  5.0, 2.0, 0.110_533_464_706_142e1],
    [  8.0, 2.0, -0.408_757_344_495_612e-1],
];
// assuming we are already in region 3
// calculate temperature given p and h
#[inline]
pub fn v_ph_3(p: Pressure, h: AvailableEnergy,) -> SpecificVolume {

    let is_region_3a = is_3a_when_in_region_3(p, h);

    if is_region_3a {
        v_ph_3a(p, h)
    } else {
        v_ph_3b(p, h)
    }

}

#[inline]
pub fn v_ph_3a(p: Pressure, h: AvailableEnergy) -> SpecificVolume {
    let v_ref = SpecificVolume::new::<cubic_meter_per_kilogram>(0.0028);
    let p_ref = Pressure::new::<megapascal>(100.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2100.0);

    let pi: f64 = (p/p_ref).into();
    let eta: f64 = (h/h_ref).into();

    // this is dimensionless volume
    let mut omega = 0.0;

    for coeffs in V_PH_SUBREGION_3A_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        omega += ni * (pi + 0.128).powi(ii as i32) * (eta - 0.727).powi(ji as i32);
    };

    return omega * v_ref;

}


/// for copying constants by hand,
/// I found it easier to use libreoffice calc
/// table 2.42 checked to be copied exactly...
const V_PH_SUBREGION_3B_COEFFS: [[f64; 3]; 30] = [
    [-12.0,0.0,-0.225_196_934_336_318E-08],
    [-12.0,1.0,0.140_674_363_313_486E-07],
    [-8.0,0.0,0.233_784_085_280_560E-05],
    [-8.0,1.0,-0.331_833_715_229_001E-04],
    [-8.0,3.0,0.107_956_778_514_318E-02],
    [-8.0,6.0,-0.271_382_067_378_863],
    [-8.0,7.0,0.107_202_262_490_333e1],
    [-8.0,8.0,-0.853_821_329_075_382],
    [-6.0,0.0,-0.215_214_194_340_526E-04],
    [-6.0,1.0,0.769_656_088_222_730E-03],
    [-6.0,2.0,-0.431_136_580_433_864E-02],
    [-6.0,5.0,0.453_342_167_309_331],
    [-6.0,6.0,-0.507_749_535_873_652],
    [-6.0,10.0,-0.100_475_154_528_389E+03],
    [-4.0,3.0,-0.219_201_924_648_793],  
    [-4.0,6.0,-0.321_087_965_668_917E+01],
    [-4.0,10.0,0.607_567_815_637_771E+03],
    [-3.0,0.0,0.557_686_450_685_932E-03],
    [-3.0,2.0,0.187_499_040_029_550],
    [-2.0,1.0,0.905_368_030_448_107E-02],
    [-2.0,2.0,0.285_417_173_048_685],
    [-1.0,0.0,0.329_924_030_996_098E-01],
    [-1.0,1.0,0.239_897_419_685_483],
    [-1.0,4.0,0.482_754_995_951_394E+01],
    [-1.0,5.0,-0.118_035_753_702_231E+02],
    [0.0,0.0,0.169_490_044_091_791],
    [1.0,0.0,-0.179_967_222_507_787E-01],
    [1.0,1.0,0.371_810_116_332_674E-01],
    [2.0,2.0,-0.536_288_335_065_096E-01],
    [2.0,6.0,0.160_697_101_092_520E+01],
];


/// eq 2.27
#[inline]
pub fn v_ph_3b(p: Pressure, h: AvailableEnergy) -> SpecificVolume {
    let v_ref = SpecificVolume::new::<cubic_meter_per_kilogram>(0.0088);
    let p_ref = Pressure::new::<megapascal>(100.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);

    let pi: f64 = (p/p_ref).into();
    let eta: f64 = (h/h_ref).into();

    // this is dimensionless volume
    let mut omega = 0.0;

    for coeffs in V_PH_SUBREGION_3B_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        omega += ni * (pi + 0.0661).powi(ii as i32) * (eta - 0.720).powi(ji as i32);
    };

    return omega * v_ref;

}
