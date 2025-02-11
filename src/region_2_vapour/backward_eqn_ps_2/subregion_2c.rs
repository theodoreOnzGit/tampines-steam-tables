
const SUBREGION_2C_BACK_COEFFS_PS: [[f64; 3]; 30] = [
    [-2.0, 0.0, 0.909_685_010_053_65e3],
    [-2.0, 1.0, 0.240_456_670_884_20e4],
    [-1.0, 0.0, -0.591_623_263_871_30e3],
    [0.0, 0.0, 0.541_454_041_280_74e3],
    [0.0, 1.0, -0.270_983_084_111_92e3],
    [0.0, 2.0, 0.979_765_250_979_26e3],
    [0.0, 3.0, -0.469_667_729_594_35e3],
    [1.0, 0.0, 0.143_992_746_047_23e2],
    [1.0, 1.0, -0.191_042_042_304_29e2],
    [1.0, 3.0, 0.532_991_671_119_71e1],
    [1.0, 4.0, -0.212_529_753_759_34e2],
    [2.0, 0.0, -0.311_473_344_137_60],
    [2.0, 1.0, 0.603_348_408_946_23],
    [2.0, 2.0, -0.427_648_397_025_09e-1],
    [3.0, 0.0, 0.581_855_972_552_59e-2],
    [3.0, 1.0, -0.145_970_082_847_53e-1],
    [3.0, 5.0, 0.566_311_756_310_27e-2],
    [4.0, 0.0, -0.761_558_645_845_77e-4],
    [4.0, 1.0, 0.224_403_429_193_32e-3],
    [4.0, 4.0, -0.125_610_950_134_13e-4],
    [5.0, 0.0, 0.633_231_326_609_34e-6],
    [5.0, 1.0, -0.205_419_896_753_75e-5],
    [5.0, 2.0, 0.364_053_703_900_82e-7],
    [6.0, 0.0, -0.297_598_977_892_15e-8],
    [6.0, 1.0, 0.101_366_185_297_63e-7],
    [7.0, 0.0, 0.599_257_196_923_51e-11],
    [7.0, 1.0, -0.206_778_701_051_64e-10],
    [7.0, 3.0, -0.208_742_781_818_86e-10],
    [7.0, 4.0, 0.101_621_668_250_89e-9],
    [7.0, 5.0, -0.164_298_282_813_47e-9],
];

use uom::si::{f64::*, pressure::megapascal, ratio::ratio, specific_heat_capacity::kilojoule_per_kilogram_kelvin, thermodynamic_temperature::kelvin};
pub(crate) fn t_ps_2c(p: Pressure, s: SpecificHeatCapacity) -> ThermodynamicTemperature{
    let p_ref = Pressure::new::<megapascal>(1.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(2.9251);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(1.0);

    let pi: f64 = (p/p_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut theta: f64 = 0.0;

    for coeffs in SUBREGION_2C_BACK_COEFFS_PS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        theta += ni * pi.powf(ii) * (2.0 - sigma).powf(ji);
    };

    return theta * t_ref;

}

