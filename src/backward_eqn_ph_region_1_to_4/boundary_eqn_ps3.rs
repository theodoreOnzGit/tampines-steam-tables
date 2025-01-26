use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal};

/// based on table 2.29
const P_S3_H_COEFFS: [[f64; 3]; 14] = [
    [0.0,0.0,6.0007364175302E-01],
    [1.0,1.0,-9.3620365484986E+00],
    [1.0,3.0,2.4659079859415E+01],
    [1.0,4.0,-1.0701422285822E+02],
    [1.0,36.0,-9.1582131580577E+13],
    [5.0,3.0,-8.6233201170066E+03],
    [7.0,0.0,-2.3583734474003E+01],
    [8.0,24.0,2.5230496938413E+17],
    [14.0,16.0,-3.8971877199772E+18],
    [20.0,16.0,-3.3377571364530E+22],
    [22.0,3.0,3.5649946963633E+10],
    [24.0,18.0,-1.4854754472064E+26],
    [28.0,8.0,3.3061151483880E+18],
    [36.0,24.0,8.1364129446783E+37],
];

#[inline]
pub fn p_s3_h(h: AvailableEnergy) -> Pressure {
    let p_ref = Pressure::new::<megapascal>(22.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2600.0);

    let eta: f64 = (h/h_ref).into();

    // this is dimensionless temperature
    let mut pi = 0.0;

    for coeffs in P_S3_H_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        pi += ni * (eta - 1.02).powi(ii as i32) * (eta - 0.608).powi(ji as i32);
    };

    return pi * p_ref;

}

