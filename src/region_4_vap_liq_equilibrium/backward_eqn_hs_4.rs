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
/// based on table 2.94
const REGION_4_BACK_COEFFS_HS: [[f64; 3]; 36] = [
    [0.0, 0.0, 0.179_882_673_606_601],
    [0.0, 3.0, -0.267_507_455_199_603],
    [0.0, 12.0, 0.116_276_722_612_600e1],
    [1.0, 0.0, 0.147_545_428_713_616],
    [1.0, 1.0, -0.512_871_635_973_248],
    [1.0, 2.0, 0.421_333_567_697_984],
    [1.0, 5.0, 0.563_749_522_189_870],
    [2.0, 0.0, 0.429_274_443_819_153],
    [2.0, 5.0, -0.335_704_552_142_140e1],
    [2.0, 8.0, 0.108_890_916_499_278e2],
    [3.0, 0.0, -0.248_483_390_456_012],
    [3.0, 2.0, 0.304_153_221_906_390],
    [3.0, 3.0, -0.494_819_763_939_905],
    [3.0, 4.0, 0.107_551_674_933_261e1],
    [4.0, 0.0, 0.733_888_415_457_688e-1],
    [4.0, 1.0, 0.140_170_545_411_085e-1],
    [5.0, 1.0, -0.106_110_975_998_808],
    [5.0, 2.0, 0.168_324_361_811_875e-1],
    [5.0, 4.0, 0.125_028_363_714_877e1],
    [5.0, 16.0, 0.101_316_840_309_509e4],
    [6.0, 6.0, -0.151_791_558_000_712e1],
    [6.0, 8.0, 0.524_277_865_990_866e2],
    [6.0, 22.0, 0.230_495_545_563_912e5],
    [8.0, 1.0, 0.249_459_806_365_456e-1],
    [10.0, 20.0, 0.210_796_467_412_137e7],
    [10.0, 36.0, 0.366_836_848_613_065e9],
    [12.0, 24.0, -0.144_814_105_365_163e9],
    [14.0, 1.0, -0.179_276_373_003_590e-2],
    [14.0, 28.0, 0.489_955_602_100_459e10],
    [16.0, 12.0, 0.471_262_212_070_518e3],
    [16.0, 32.0, -0.829_294_390_198_652e11],
    [18.0, 14.0, -0.171_545_662_263_191e4],
    [18.0, 22.0, 0.355_777_682_973_575e7],
    [18.0, 36.0, 0.586_062_760_258_436e12],
    [20.0, 24.0, -0.129_887_635_078_195e8],
    [28.0, 36.0, 0.317_247_449_371_057e11],
];
