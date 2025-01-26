use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, ratio::ratio};
/// based on eq 2.25
#[inline]
pub fn h_3a3b_backwards_ph_boundary(p: Pressure) -> AvailableEnergy {
    let p_ref = Pressure::new::<megapascal>(1.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1.0);

    let pi: f64 = (p/p_ref).get::<ratio>();

    let n1: f64 = 0.201_464_004_206_875e4;
    let n2: f64 = 0.374_696_550_136_983e1;
    let n3: f64 = -0.219_921_901_054_187_e-1;
    let n4: f64 = 0.875_131_686_009_950e-4;

    let eta = n1 
        + n2 * pi
        + n3 * pi.powi(2)
        + n4 * pi.powi(3);

    return eta * h_ref;

}
