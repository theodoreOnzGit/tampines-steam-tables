use crate::region_1_subcooled_liquid::REGION_1_BACK_COEFFS_PS;

pub(crate) fn t_ps_1_kelvin(p: f64, s: f64) -> f64 {
    let sig = sigma_1_back_ps(s);
    let pi = pi_1_back_ps(p);
    let mut sum = 0.0;
    for coefficient in REGION_1_BACK_COEFFS_PS {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += ni * pi.powi(ii) * (sig + 2.0).powi(ji);
    }
    sum
}


/// Returns the region-1 pi for backwards calculations
/// Pressure is assumed to be in Pa
pub(crate) fn pi_1_back_ps(p_pascals: f64) -> f64 {
    p_pascals / 1e6
}

/// Returns the region-1 sigma for backwards calculations
/// Entropy is assumed to be in kJ/kg.K
pub(crate) fn sigma_1_back_ps(s_kj_per_kg_k: f64) -> f64 {
    s_kj_per_kg_k
}
