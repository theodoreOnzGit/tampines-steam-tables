use uom::si::{f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin};


/// based on table 2.29
const P_S3_S_COEFFS: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.639_767_553_612_785],
    [1.0, 1.0, -0.129_727_445_396_014e2],
    [1.0, 32.0, -0.224_595_125_848_403e16],
    [4.0, 7.0, 0.177_466_741_801_846e7],
    [12.0, 4.0, 0.717_079_349_571_538e10],
    [12.0, 14.0, -0.378_829_107_169_011e18],
    [16.0, 36.0, -0.955_586_736_431_328e35],
    [24.0, 10.0, 0.187_269_814_676_188e24],
    [28.0, 0.0, 0.119_254_746_466_473e12],
    [32.0, 18.0, 0.110_649_277_244_882e37],
];

#[inline]
pub fn p_s3_s(s: SpecificHeatCapacity) -> Pressure {
    let p_ref = Pressure::new::<megapascal>(22.0);
    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.2);

    let sigma: f64 = (s/s_ref).into();

    // this is dimensionless temperature
    let mut pi = 0.0;

    for coeffs in P_S3_S_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (sigma - 1.03).powi(ii as i32) * (sigma - 0.699).powi(ji as i32);
    };

    return pi * p_ref;

}

