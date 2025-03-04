use uom::si::pressure::megapascal;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

const REGION_1_BACK_COEFFS_HS: [[f64; 3]; 19] = [
    [0.0, 0.0, -0.691_997_014_660_582],
    [0.0, 1.0, -0.183_612_548_787_560e2],
    [0.0, 2.0, -0.928_332_409_297_335e1],
    [0.0, 4.0, 0.659_639_569_909_906e2],
    [0.0, 5.0, -0.162_060_388_912_024e2],
    [0.0, 6.0, 0.450_620_017_338_667e3],
    [0.0, 8.0, 0.854_680_678_224_170e3],
    [0.0, 14.0, 0.607_523_214_001_162e4],
    [1.0, 0.0, 0.326_487_682_621_856e2],
    [1.0, 1.0, -0.269_408_844_582_931e2],
    [1.0, 4.0, -0.319_947_848_334_300e3],
    [1.0, 6.0, -0.928_354_307_043_320e3],
    [2.0, 0.0, 0.303_634_537_455_249e2],
    [2.0, 1.0, -0.650_540_422_444_146e2],
    [2.0, 10.0, -0.430_991_316_516_130e4],
    [3.0, 4.0, -0.747_512_324_096_068e3],
    [4.0, 1.0, 0.730_000_345_529_245e3],
    [4.0, 4.0, 0.114_284_032_269_021e4],
    [5.0, 0.0, -0.436_407_041_874_559e3],
];

pub fn p_hs_1(h: AvailableEnergy, s: SpecificHeatCapacity) -> Pressure {
    let p_ref = Pressure::new::<megapascal>(100.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(3400.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(7.6);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut pi: f64 = 0.0;


    for coefficient in REGION_1_BACK_COEFFS_HS {
        let ii = coefficient[0];
        let ji = coefficient[1];
        let ni = coefficient[2];
        pi  += ni * (eta + 0.05).powf(ii) * (sigma + 0.05).powf(ji);
    }

    return pi * p_ref;
}
