use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::thermodynamic_temperature::kelvin;

// assuming we are already in region 3
// calculate temperature given p and h
#[inline]
pub fn tsat_hs_4(h: AvailableEnergy, s: SpecificHeatCapacity) -> ThermodynamicTemperature {
    let t_ref = ThermodynamicTemperature::new::<kelvin>(550.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.2);

    let eta: f64 = (h/h_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut theta: f64 = 0.0;

    for coeffs in REGION_4_BACK_COEFFS_HS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        theta += ni * (eta - 0.119).powf(ii) * (sigma - 1.07).powf(ji);
    };

    return theta * t_ref;

}
/// based on table 2.88
const REGION_4_BACK_COEFFS_HS: [[f64; 3]; 36] = [
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
    [32.0, 28.0, 0.377_121_605_943_324e41],
    [32.0, 28.0, 0.377_121_605_943_324e41],
    [32.0, 28.0, 0.377_121_605_943_324e41],
];
