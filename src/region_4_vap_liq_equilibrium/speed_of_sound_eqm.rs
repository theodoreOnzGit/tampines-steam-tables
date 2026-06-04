use uom::si::f64::*;
use uom::si::pressure::megapascal;
use uom::si::ratio::ratio;
use uom::si::temperature_interval::degree_celsius;
use uom::si::thermodynamic_temperature::kelvin;

use crate::interfaces::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::region_1_subcooled_liquid::{
    v_tp_1,
    h_tp_1,
    s_tp_1,
};
use crate::region_2_vapour::{
    v_tp_2,
    h_tp_2,
    s_tp_2,
};
use crate::region_4_vap_liq_equilibrium::
    sat_temp_4
;

/// this is a more complicated version that makes use of derivatives 
/// directly based on thermodynamic calculus 
///
/// AI generated...
///
pub fn w_ps_eqm_region4_kieffer(p: Pressure, s: SpecificHeatCapacity) -> Velocity {
    
    //let r = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(461.526);
    let t_sat = sat_temp_4(p);
    
    // --- saturation properties ---
    let v_f = v_tp_1(t_sat, p);  // specific volume liquid
    let v_g = v_tp_2(t_sat, p);  // specific volume vapour
    let s_f = s_tp_1(t_sat, p);  // entropy liquid
    let s_g = s_tp_2(t_sat, p);  // entropy vapour
    let h_f = h_tp_1(t_sat, p);  // enthalpy liquid
    let h_g = h_tp_2(t_sat, p);  // enthalpy vapour
    let l   = h_g - h_f;         // latent heat
    
    // quality from entropy
    let x: f64 = ((s - s_f) / (s_g - s_f)).get::<ratio>();
    let eta = x; // Kieffer uses eta for mass fraction of steam
    
    // mixture specific volume
    let v_mix = v_f * (1.0 - eta) + v_g * eta;
    
    // --- dT_sat/dp via finite difference ---
    let dp = p * 1e-6_f64;
    let dt_kelvin = sat_temp_4(p + dp).get::<kelvin>() - sat_temp_4(p - dp).get::<kelvin>();
    let _dt_dp = TemperatureInterval::new::<degree_celsius>(dt_kelvin) / (2.0 * dp);
    
    // --- dv_f/dp|_sat via finite difference along saturation curve ---
    let v_f_plus  = v_tp_1(sat_temp_4(p + dp), p + dp);
    let v_f_minus = v_tp_1(sat_temp_4(p - dp), p - dp);
    let dv_f_dp = (v_f_plus - v_f_minus) / (2.0 * dp);
    
    // --- dv_g/dp|_sat via finite difference along saturation curve ---
    let v_g_plus  = v_tp_2(sat_temp_4(p + dp), p + dp);
    let v_g_minus = v_tp_2(sat_temp_4(p - dp), p - dp);
    let dv_g_dp = (v_g_plus - v_g_minus) / (2.0 * dp);
    
    // --- dh_f/dp|_sat via finite difference ---
    let h_f_plus  = h_tp_1(sat_temp_4(p + dp), p + dp);
    let h_f_minus = h_tp_1(sat_temp_4(p - dp), p - dp);
    let dh_f_dp = (h_f_plus - h_f_minus) / (2.0 * dp);
    
    // --- dL/dp|_sat via finite difference ---
    let l_plus  = h_tp_2(sat_temp_4(p + dp), p + dp) 
                - h_tp_1(sat_temp_4(p + dp), p + dp);
    let l_minus = h_tp_2(sat_temp_4(p - dp), p - dp) 
                - h_tp_1(sat_temp_4(p - dp), p - dp);
    let dl_dp = (l_plus - l_minus) / (2.0 * dp);
    
    // --- Kieffer eq. 28 ---
    // -c²/V² = (1-η)(dV_f/dP)_sat + η(dV_g/dP)_sat 
    //        + (V_g - V_f) * [V/L - (1/L)(dH_f/dP)_sat - (η/L)(dL/dP)]
    
    let term1 = (1.0 - eta) * dv_f_dp;
    let term2 = eta * dv_g_dp;
    let bracket = v_mix / l 
                - dh_f_dp / l 
                - eta * dl_dp / l;
    let term3 = (v_g - v_f) * bracket;
    
    let neg_c_sq_over_v_sq = term1 + term2 + term3;
    
    // c² = -V² / (neg_c_sq_over_v_sq / (-1))
    // since neg_c_sq_over_v_sq should be negative:
    let c_sq = v_mix * v_mix / (-neg_c_sq_over_v_sq);
    
    c_sq.sqrt()
}

/// this is a simpler version that makes use of derivatives 
/// that makes use of derivatives
/// AI generated
pub fn w_ps_eqm_region4_finite_diff_vol(p: Pressure, s: SpecificHeatCapacity) -> Velocity {
    
    // guard against going below minimum steam table pressure
    let p_min = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let dp = p * 1e-4_f64;
    let p_plus  = p + dp;
    let p_minus = if p - dp > p_min { p - dp } else { p_min };

    let v       = v_ps_eqm(p,       s);
    let v_plus  = v_ps_eqm(p_plus,  s);
    let v_minus = v_ps_eqm(p_minus, s);

    let dp_actual = p_plus - p_minus;
    let dv_dp_s = (v_plus - v_minus) / dp_actual;

    // dv_dp_s should be negative in two-phase region
    // c = v * sqrt(1 / (-dv/dp|_s))
    let c_hem: Velocity = v * (dv_dp_s * -1.0).recip().sqrt();

    c_hem
}

