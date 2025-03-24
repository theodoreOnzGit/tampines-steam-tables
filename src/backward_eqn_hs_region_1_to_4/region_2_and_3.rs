
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::thermodynamic_temperature::kelvin;

/// this is for eq 2.44 on page 84
/// based on table 2.73
const TB23_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 25] = [
    [-12.0, 10.0, 0.629_096_260_829_810e-3],
    [-10.0, 8.0, -0.823_453_502_583_165e-3],
    [-8.0, 3.0, 0.515_446_951_519_474e-7],
    [-4.0, 4.0, -0.117_565_945_784_945e1],
    [-3.0, 3.0, 0.348_519_684_726_192e1],
    [-2.0, -6.0, -0.507_837_382_408_313e-11],
    [-2.0, 2.0, -0.284_637_670_005_479e1],
    [-2.0, 3.0, -0.236_092_263_939_673e1],
    [-2.0, 4.0, 0.601_492_324_973_779e1],
    [0.0, 0.0, 0.148_039_650_824_546e1],
    [1.0, -3.0, 0.360_075_182_221_907e-3],
    [1.0, -2.0, -0.126_700_045_009_952e-1],
    [1.0, 10.0, -0.122_184_332_521_413e7],
    [3.0, -2.0, 0.149_276_502_463_272],
    [3.0, -1.0, 0.698_733_471_798_484],
    [5.0, -5.0, -0.252_207_040_114_321e-1],
    [6.0, -6.0, 0.147_151_930_985_213e-1],
    [6.0, -3.0, -0.108_618_917_681_849e1],
    [8.0, -8.0, -0.936_875_039_816_322e-3],
    [8.0, -2.0, 0.819_877_897_570_217e2],
    [8.0, -1.0, -0.182_041_861_521_835e3],
    [12.0, -12.0, 0.261_907_376_402_688e-5],
    [12.0, -1.0, -0.291_626_417_025_961e5],
    [14.0, -12.0, 0.140_660_774_926_165e-4],
    [14.0, 1.0, 0.783_237_062_349_385e7],
];


/// this function represents the saturated liquid line
/// for hs flashing between region 1 and region 4
pub fn tb23_s_boundary_enthalpy(
    s: SpecificHeatCapacity,
    h: AvailableEnergy) -> ThermodynamicTemperature {

    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.3);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(3000.0);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(900.0);
    let sigma: f64 = (s/s_ref).get::<ratio>();
    let eta: f64 = (h/h_ref).get::<ratio>();

    let mut theta: f64 = 0.0;

    for coeffs in TB23_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        theta += ni * (eta - 0.727).powf(ii) * (sigma - 0.864).powf(ji);
    }

    return t_ref * theta;

}



