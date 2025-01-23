const REGION_1_BACK_COEFFS_PH: [[f64; 3]; 20] = [
    [0.0, 0.0, -0.23872489924521e+3],
    [0.0, 1.0, 0.40421188637945e+3],
    [0.0, 2.0, 0.11349746881718e+3],
    [0.0, 6.0, -0.58457616048039e+1],
    [0.0, 22.0, -0.15285482413140e-3],
    [0.0, 32.0, -0.10866707695377e-5],
    [1.0, 0.0, -0.13391744872602e+2],
    [1.0, 1.0, 0.43211039183559e+2],
    [1.0, 2.0, -0.54010067170506e+2],
    [1.0, 3.0, 0.30535892203916e+2],
    [1.0, 4.0, -0.65964749423638e+1],
    [1.0, 10.0, 0.93965400878363e-2],
    [1.0, 32.0, 0.11573647505340e-6],
    [2.0, 10.0, -0.25858641282073e-4],
    [2.0, 32.0, -0.40644363084799e-8],
    [3.0, 10.0, 0.66456186191635e-7],
    [3.0, 32.0, 0.80670734103027e-10],
    [4.0, 32.0, -0.93477771213947e-12],
    [5.0, 32.0, 0.58265442020601e-14],
    [6.0, 32.0, -0.15020185953503e-16],
];

use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, thermodynamic_temperature::degree_celsius};
/// Returns the region-1 eta for backwards calculations
/// Enthalpy is assumed to be in kJ/kg
pub fn eta_1_back(h: AvailableEnergy) -> f64 {
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2500.0);
    (h / h_ref).into()
}

/// Returns the region-1 pi for backwards calculations
pub fn pi_1_back(p: Pressure) -> f64 {
    let p_ref = Pressure::new::<megapascal>(1.0);
    (p / p_ref).into()
}

/// Returns the region-1 backward correlation for T(p,h)
///
/// the reference temperature is 1K
pub fn t_ph_1(p: Pressure, h: AvailableEnergy,) -> ThermodynamicTemperature {
    let t = ThermodynamicTemperature::new::<degree_celsius>(
        1.0 * theta_ph_1(p, h)
    );

    t

}

/// Returns the region-1 backward correlation for theta = T/T* (p,h)
#[inline]
pub fn theta_ph_1(p: Pressure, h: AvailableEnergy) -> f64 {
    let eta = eta_1_back(h);
    let pi = pi_1_back(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_BACK_COEFFS_PH {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii) * (eta + 1.0).powi(ji);
    }
    sum
}

