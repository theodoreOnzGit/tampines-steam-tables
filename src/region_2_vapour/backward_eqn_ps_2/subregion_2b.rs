
const SUBREGION_2B_BACK_COEFFS_PS: [[f64; 3]; 44] = [
    [-6.0, 0.0, 0.316_876_650_834_97e6],
    [-6.0, 11.0, 0.208_641_758_818_58e2],
    [-5.0, 0.0, -0.398_593_998_035_99e6],
    [-5.0, 11.0, -0.218_160_585_188_77e2],
    [-4.0, 0.0, 0.223_697_851_942_42e6],
    [-4.0, 1.0, -0.278_417_034_458_17e4],
    [-4.0, 11.0, 0.992_074_360_714_80e1],
    [-3.0, 0.0, -0.751_975_122_991_57e5],
    [-3.0, 1.0, 0.297_086_059_511_58e4],
    [-3.0, 11.0, -0.344_068_785_485_26e1],
    [-3.0, 12.0, 0.388_155_642_491_15],
    [-2.0, 0.0, 0.175_112_950_857_50e5],
    [-2.0, 1.0, -0.142_371_128_544_49e4],
    [-2.0, 6.0, 0.109_438_033_641_67e1],
    [-2.0, 10.0, 0.899_716_193_084_95],
    [-1.0, 0.0, -0.337_597_400_989_58e4],
    [-1.0, 1.0, 0.471_628_858_183_55e3],
    [-1.0, 5.0, -0.191_882_419_936_79e1],
    [-1.0, 8.0, 0.410_785_804_921_96],
    [-1.0, 9.0, -0.334_653_781_720_97],
    [0.0, 0.0, 0.138_700_347_775_05e4],
    [0.0, 1.0, -0.406_633_261_952_38e3],
    [0.0, 2.0, 0.417_273_471_596_10e2],
    [0.0, 4.0, 0.219_325_494_345_32e1],
    [0.0, 5.0, -0.103_200_500_090_77e1],
    [0.0, 6.0, 0.358_829_435_167_03],
    [0.0, 9.0, 0.525_114_537_260_66e-2],
    [1.0, 0.0, 0.128_389_164_507_05e2],
    [1.0, 1.0, -0.286_424_372_193_81e1],
    [1.0, 2.0, 0.569_126_836_648_55],
    [1.0, 3.0, -0.999_629_545_849_31e-1],
    [1.0, 7.0, -0.326_320_377_784_59e-2],
    [1.0, 8.0, 0.233_209_225_767_23e-3],
    [2.0, 0.0, -0.153_348_098_574_50],
    [2.0, 1.0, 0.290_722_882_399_02e-1],
    [2.0, 5.0, 0.375_347_027_411_67e-3],
    [3.0, 0.0, 0.172_966_917_024_11e-2],
    [3.0, 1.0, -0.385_560_508_445_04e-3],
    [3.0, 3.0, -0.350_177_122_926_08e-4],
    [4.0, 0.0, -0.145_663_936_314_92e-4],
    [4.0, 1.0, 0.564_208_572_672_69e-5],
    [5.0, 0.0, 0.412_861_500_746_05e-7],
    [5.0, 1.0, -0.206_846_711_188_24e-7],
    [5.0, 2.0, 0.164_093_936_747_25e-8],
];

use uom::si::{f64::*, pressure::megapascal, ratio::ratio, specific_heat_capacity::kilojoule_per_kilogram_kelvin, thermodynamic_temperature::kelvin};
pub(crate) fn t_ps_2b(p: Pressure, s: SpecificHeatCapacity) -> ThermodynamicTemperature{
    let p_ref = Pressure::new::<megapascal>(1.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(0.7853);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(1.0);

    let pi: f64 = (p/p_ref).get::<ratio>();
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut theta: f64 = 0.0;

    for coeffs in SUBREGION_2B_BACK_COEFFS_PS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        theta += ni * pi.powf(ii) * (10.0 - sigma).powf(ji);
    };

    return theta * t_ref;

}

