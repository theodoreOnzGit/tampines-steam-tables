
/// Returns the region-2 tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn tau_2(t: f64) -> f64 {
    540.0 / t
}

/// Returns the region-2 pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn pi_2(p: f64) -> f64 {
    p / 1e6
}
