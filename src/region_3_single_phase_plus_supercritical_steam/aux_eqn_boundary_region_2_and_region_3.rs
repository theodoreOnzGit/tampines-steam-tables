use uom::si::{f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

/// boundary equation between region 2 and 3
/// note that points ON this line belong to region 2
pub fn p_boundary_2_3(t: ThermodynamicTemperature) -> Pressure {
    let n: [f64; 5] = [
        0.34805185628969e3,
        -0.11671859879975e1,
        0.10192970039326e-2,
        0.57254459862746e3,
        0.13918839778870e2,
    ];
    let p_ref = Pressure::new::<megapascal>(1.0);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(1.0);
    // theta is dimensionless temp
    let theta: f64 = (t/t_ref).into();
    let dimensionless_pressure = n[0] + n[1] * theta + n[2] * theta.powi(2);
    
    return p_ref * dimensionless_pressure;

}
/// boundary equation between region 2 and 3
/// note that points ON this line belong to region 2
pub fn t_boundary_2_3(p: Pressure) -> ThermodynamicTemperature {
    let n: [f64; 5] = [
        0.34805185628969e3,
        -0.11671859879975e1,
        0.10192970039326e-2,
        0.57254459862746e3,
        0.13918839778870e2,
    ];
    let p_ref = Pressure::new::<megapascal>(1.0);
    let t_ref = ThermodynamicTemperature::new::<kelvin>(1.0);
    let dimensionless_p: f64 = (p/p_ref).into();
    // theta is dimensionless temp
    let theta = n[3] + (((dimensionless_p) - n[4]) / n[2]).sqrt();
    return t_ref * theta;
    
}
