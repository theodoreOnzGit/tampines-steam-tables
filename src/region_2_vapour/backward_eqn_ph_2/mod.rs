use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::{megapascal, pascal}, thermodynamic_temperature::kelvin};
#[inline]
pub fn t_ph_2(p: Pressure, h: AvailableEnergy) -> ThermodynamicTemperature {

    let p_ref = Pressure::new::<megapascal>(1.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2000.0);
    let pi: f64 = (p / p_ref).into();
    let eta: f64 = (h / h_ref).into();
    let p_2b2c = p_2b2c(h);

    match p.get::<pascal>() {
        pres if (0.0..=4.0e6).contains(&pres) => t_ph_2a(pi, eta),
        pres if (4.0e6..=100.0e6).contains(&pres) && pres < p_2b2c.get::<pascal>() => t_ph_2b(pi, eta),
        _ => t_ph_2c(pi, eta),
    }
}

/// eqn for determining pressure boundary between subregion 2b and 2c
/// using dimensionless enthalpy eta
#[inline]
pub fn p_2b2c(h: AvailableEnergy) -> Pressure {

    let p_ref = Pressure::new::<megapascal>(1.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0);
    let eta: f64 = (h/h_ref).into();
    let n1 = 0.90584278514723e3;
    let n2 = -0.67955786399241;
    let n3 = 0.12809002730136e-3;
    let p_2b2c = (n1 +  n2 * eta +  n3 * eta.powi(2)) * p_ref;

    p_2b2c
}

/// eqn for determining enthalpy boundary between subregion 2b and 2c
/// using dimensionless pressure pi
#[inline] 
pub fn h_2b2c(p: Pressure) -> AvailableEnergy {
    let p_ref = Pressure::new::<megapascal>(1.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0);
    let pi: f64 = (p/p_ref).into();
    let n3 = 0.12809002730136e-3;
    let n4 = 0.265_265_719_084_28e4;
    let n5 = 0.452_575_789_059_48e1;

    let h_2b2c = n4 + ( (pi-n5)/n3 ).sqrt();

    return h_2b2c * h_ref;

}

#[inline]
pub fn t_ph_2a(pi: f64, eta: f64) -> ThermodynamicTemperature {
    let i: [i32; 34] = [
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5,
        5, 6, 6, 7,
    ];
    let j: [i32; 34] = [
        0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12,
        32, 44, 32, 36, 42, 34, 44, 28,
    ];
    let n: [f64; 34] = [
        0.10898952318288e4,
        0.84951654495535e3,
        -0.10781748091826e3,
        0.33153654801263e2,
        -0.74232016790248e1,
        0.11765048724356e2,
        0.18445749355790e1,
        -0.41792700549624e1,
        0.62478196935812e1,
        -0.17344563108114e2,
        -0.20058176862096e3,
        0.27196065473796e3,
        -0.45511318285818e3,
        0.30919688604755e4,
        0.25226640357872e6,
        -0.61707422868339e-2,
        -0.31078046629583,
        0.11670873077107e2,
        0.12812798404046e9,
        -0.98554909623276e9,
        0.28224546973002e10,
        -0.35948971410703e10,
        0.17227349913197e10,
        -0.13551334240775e5,
        0.12848734664650e8,
        0.13865724283226e1,
        0.23598832556514e6,
        -0.13105236545054e8,
        0.73999835474766e4,
        -0.55196697030060e6,
        0.37154085996233e7,
        0.19127729239660e5,
        -0.41535164835634e6,
        -0.62459855192507e2,
    ];

    // Calculate T
    let x: [usize; 34] = core::array::from_fn(|i| i + 1);
    let theta: f64 = x.into_iter()
        .map(|x| n[x - 1] * pi.powi(i[x - 1]) * (eta - 2.1).powi(j[x - 1]))
        .sum();

    let t_ref_kelvin: f64 = 1.0;
    return ThermodynamicTemperature::new::<kelvin>(theta * t_ref_kelvin);
}

#[inline]
pub fn t_ph_2b(pi: f64, eta: f64) -> ThermodynamicTemperature {
    let i: [i32; 38] = [
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 6, 7, 7, 9, 9,
    ];
    let j: [i32; 38] = [
        0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2,
        12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40,
    ];
    let n: [f64; 38] = [
        0.14895041079516e4,
        0.74307798314034e3,
        -0.97708318797837e2,
        0.24742464705674e1,
        -0.63281320016026,
        0.11385952129658e1,
        -0.47811863648625,
        0.85208123431544e-2,
        0.93747147377932,
        0.33593118604916e1,
        0.33809355601454e1,
        0.16844539671904,
        0.73875745236695,
        -0.47128737436186,
        0.15020273139707,
        -0.21764114219750e-2,
        -0.21810755324761e-1,
        -0.10829784403677,
        -0.46333324635812e-1,
        0.71280351959551e-4,
        0.11032831789999e-3,
        0.18955248387902e-3,
        0.30891541160537e-2,
        0.13555504554949e-2,
        0.28640237477456e-6,
        -0.10779857357512e-4,
        -0.76462712454814e-4,
        0.14052392818316e-4,
        -0.31083814331434e-4,
        -0.10302738212103e-5,
        0.28217281635040e-6,
        0.12704902271945e-5,
        0.73803353468292e-7,
        -0.11030139238909e-7,
        -0.81456365207833e-13,
        -0.25180545682962e-10,
        -0.17565233969407e-17,
        0.86934156344163e-14,
    ];

    // Calculate T
    let x: [usize; 38] = core::array::from_fn(|i| i + 1);

    let theta: f64 = x.into_iter()
        .map(|x| n[x - 1] * (pi - 2.0).powi(i[x - 1]) * (eta - 2.6).powi(j[x - 1]))
        .sum();
    let t_ref_kelvin: f64 = 1.0;
    return ThermodynamicTemperature::new::<kelvin>(theta * t_ref_kelvin);
}

#[inline]
pub fn t_ph_2c(pi: f64, eta: f64) -> ThermodynamicTemperature {
    let i: [i32; 23] = [
        -7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6,
    ];
    let j: [i32; 23] = [
        0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22,
    ];
    let n: [f64; 23] = [
        -0.32368398555242e13,
        0.73263350902181e13,
        0.35825089945447e12,
        -0.58340131851590e12,
        -0.10783068217470e11,
        0.20825544563171e11,
        0.61074783564516e6,
        0.85977722535580e6,
        -0.25745723604170e5,
        0.31081088422714e5,
        0.12082315865936e4,
        0.48219755109255e3,
        0.37966001272486e1,
        -0.10842984880077e2,
        -0.45364172676660e-1,
        0.14559115658698e-12,
        0.11261597407230e-11,
        -0.17804982240686e-10,
        0.12324579690832e-6,
        -0.11606921130984e-5,
        0.27846367088554e-4,
        -0.59270038474176e-3,
        0.12918582991878e-2,
    ];

    // Calculate T
    let x: [usize; 23] = core::array::from_fn(|i| i + 1);
    let theta: f64 = x.into_iter()
        .map(|x| n[x - 1] * (pi + 25.0).powi(i[x - 1]) * (eta - 1.8).powi(j[x - 1]))
        .sum();
    let t_ref_kelvin: f64 = 1.0;
    return ThermodynamicTemperature::new::<kelvin>(theta * t_ref_kelvin);
}

