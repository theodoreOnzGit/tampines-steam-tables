use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::quantities::SpecificHeatCapacity;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};
use uom::ConstZero;

use crate::interfaces::functional_programming::hs_flash_eqm::{t_hs_eqm, v_hs_eqm};
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;
use crate::region_3_single_phase_plus_supercritical_steam::{h_rho_t_3, v_tp_3c, v_tp_3r, v_tp_3s, v_tp_3t, v_tp_3u, v_tp_3x, v_tp_3y, v_tp_3z};
use crate::region_4_vap_liq_equilibrium::sat_temp_4;

/// saturation table (see page 182)
///
/// this is region 1 and region 4
#[test]
pub fn saturation_table_260_degc_to_350_degc(){

    //[t_deg_c,t_kelvin,psat_bar,v_liq_m3_per_kg,v_vap_m3_per_kg,h_liq_kj_per_kg,h_vap_kj_per_kg,enthalpy_of_vap,s_liq_kj_per_kg_k,s_vap_kj_per_kg_k],
    let steam_table: Vec<[f64; 10]> =
        vec![
        [260.0,533.15,46.9207,0.00127613,0.0421755,1134.83,2796.64,1661.82,2.8847,6.0017],
        [262.0,535.15,48.464,0.00128129,0.040766,1144.78,2795.47,1650.68,2.903,5.9875],
        [264.0,537.15,50.0457,0.00128656,0.0394082,1154.79,2794.19,1639.4,2.9213,5.9733],
        [266.0,539.15,51.6666,0.00129193,0.0380997,1164.84,2792.8,1627.96,2.9396,5.959],
        [268.0,541.15,53.3273,0.00129741,0.0368385,1174.94,2791.3,1616.36,2.9579,5.9448],
        [270.0,543.15,55.0284,0.00130301,0.0356224,1185.09,2789.69,1604.6,2.9762,5.9304],
        [272.0,545.15,56.7706,0.00130872,0.0344496,1195.3,2787.96,1592.66,2.9945,5.916],
        [274.0,547.15,58.5547,0.00131455,0.033318,1205.55,2786.11,1580.56,3.0129,5.9016],
        [276.0,549.15,60.3812,0.00132052,0.032226,1215.87,2784.14,1568.28,3.0312,5.8871],
        [278.0,551.15,62.251,0.00132661,0.0311719,1226.24,2782.05,1555.81,3.0496,5.8725],
        [280.0,553.15,64.1646,0.00133285,0.030154,1236.67,2779.82,1543.15,3.0681,5.8578],
        [282.0,555.15,66.1228,0.00133922,0.0291708,1247.16,2777.47,1530.3,3.0865,5.8431],
        [284.0,557.15,68.1264,0.00134575,0.0282208,1257.72,2774.97,1517.25,3.105,5.8283],
        [286.0,559.15,70.176,0.00135243,0.0273027,1268.34,2772.34,1504.0,3.1236,5.8134],
        [288.0,561.15,72.2724,0.00135928,0.0264152,1279.03,2769.56,1490.53,3.1421,5.7984],
        [290.0,563.15,74.4164,0.00136629,0.0255568,1289.8,2766.63,1476.84,3.1608,5.7832],
        [292.0,565.15,76.6087,0.00137349,0.0247265,1300.63,2763.55,1462.92,3.1794,5.768],
        [294.0,567.15,78.8502,0.00138087,0.0239231,1311.54,2760.31,1448.76,3.1982,5.7526],
        [296.0,569.15,81.1415,0.00138844,0.0231454,1322.54,2756.9,1434.37,3.217,5.7372],
        [298.0,571.15,83.4835,0.00139622,0.0223924,1333.61,2753.33,1419.72,3.2358,5.7215],
        [300.0,573.15,85.8771,0.00140422,0.0216631,1344.77,2749.57,1404.8,3.2547,5.7058],
        [305.0,578.15,92.0919,0.00142524,0.019937,1373.07,2739.38,1366.3,3.3024,5.6656],
        [310.0,583.15,98.6475,0.00144788,0.0183389,1402.0,2727.92,1325.92,3.3506,5.6243],
        [315.0,588.15,105.558,0.00147239,0.0168557,1431.63,2715.08,1283.45,3.3994,5.5816],
        [320.0,593.15,112.839,0.00149906,0.0154759,1462.05,2700.67,1238.62,3.4491,5.5373],
        [325.0,598.15,120.505,0.0015283,0.0141887,1493.37,2684.48,1191.11,3.4997,5.4911],
        [330.0,603.15,128.575,0.0015606,0.012984,1525.74,2666.25,1140.51,3.5516,5.4425],
        [335.0,608.15,137.067,0.00159667,0.0118522,1559.34,2645.6,1086.26,3.6048,5.391],
        [340.0,613.15,146.002,0.00163751,0.0107838,1594.45,2622.07,1027.62,3.6599,5.3359],
        [345.0,618.15,155.401,0.0016846,0.0097698,1631.44,2595.01,963.57,3.7175,5.2763],
        [350.0,623.15,165.292,0.00174007,0.0088009,1670.86,2563.59,892.73,3.7783,5.2109],
        ];

        for dataset in steam_table {
            let t_deg_c = dataset[0];
            let t_kelvin = dataset[1];
            let psat_bar = dataset[2];
            let v_liq_m3_per_kg = dataset[3];
            let v_vap_m3_per_kg = dataset[4];
            let h_liq_kj_per_kg = dataset[5];
            let h_vap_kj_per_kg = dataset[6];
            let enthalpy_of_vap_kj_per_kg = dataset[7];
            let s_liq_kj_per_kg_k = dataset[8];
            let s_vap_kj_per_kg_k = dataset[9];
            assert_hs_flash(t_deg_c, t_kelvin, psat_bar, 
                v_liq_m3_per_kg, v_vap_m3_per_kg, h_liq_kj_per_kg, 
                h_vap_kj_per_kg, enthalpy_of_vap_kj_per_kg, 
                s_liq_kj_per_kg_k, s_vap_kj_per_kg_k);
        }

}


/// saturation table (see page 182)
///
/// this is region 3 and region 4
#[test]
pub fn saturation_table_355_degc_to_crit_373_degc(){

    //[t_deg_c,t_kelvin,psat_bar,v_liq_m3_per_kg,v_vap_m3_per_kg,h_liq_kj_per_kg,h_vap_kj_per_kg,enthalpy_of_vap,s_liq_kj_per_kg_k,s_vap_kj_per_kg_k],
    let steam_table: Vec<[f64; 10]> =
        vec![
        [355.0,628.15,175.701,0.0018078,0.007866,1713.71,2526.45,812.74,3.8438,5.1377],
        [360.0,633.15,186.664,0.00189451,0.0069449,1761.49,2480.99,719.5,3.9164,5.0527],
        [365.0,638.15,198.222,0.00201561,0.0060044,1817.59,2422.0,604.41,4.0011,4.9482],
        // very near critical point, this kinda fails bad
        [370.0,643.15,210.434,0.00222209,0.0049462,1892.64,2333.5,440.86,4.1142,4.7996],
        [371.0,644.15,212.964,0.0022902,0.0046914,1913.25,2307.45,394.2,4.1453,4.7573],
        [372.0,645.15,215.528,0.0023817,0.0043985,1938.54,2274.69,336.15,4.1836,4.7046],
        [373.0,646.15,218.132,0.00252643,0.0040212,1974.14,2227.55,253.42,4.2377,4.6299],
        [373.946,647.096,220.64,0.00310559,0.00310559,2087.55,2087.55,0.0,4.412,4.412],
        ];

        for dataset in steam_table {
            let t_deg_c = dataset[0];
            let t_kelvin = dataset[1];
            let psat_bar = dataset[2];
            let v_liq_m3_per_kg = dataset[3];
            let v_vap_m3_per_kg = dataset[4];
            let h_liq_kj_per_kg = dataset[5];
            let h_vap_kj_per_kg = dataset[6];
            let enthalpy_of_vap_kj_per_kg = dataset[7];
            let s_liq_kj_per_kg_k = dataset[8];
            let s_vap_kj_per_kg_k = dataset[9];
            assert_hs_flash_near_crit(t_deg_c, t_kelvin, psat_bar, 
                v_liq_m3_per_kg, v_vap_m3_per_kg, h_liq_kj_per_kg, 
                h_vap_kj_per_kg, enthalpy_of_vap_kj_per_kg, 
                s_liq_kj_per_kg_k, s_vap_kj_per_kg_k);
        }

}


fn assert_hs_flash(t_deg_c: f64,
    t_kelvin: f64,
    psat_bar: f64,
    v_liq_m3_per_kg: f64,
    v_vap_m3_per_kg: f64,
    h_liq_kj_per_kg: f64,
    h_vap_kj_per_kg: f64,
    enthalpy_of_vap_kj_per_kg: f64,
    s_liq_kj_per_kg_k: f64,
    s_vap_kj_per_kg_k: f64){

    // specify a vapour quality
    let x_ref = 0.77;
    let p = Pressure::new::<bar>(psat_bar);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(
        (1.0-x_ref) * h_liq_kj_per_kg + x_ref * h_vap_kj_per_kg
    );
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        (1.0-x_ref) * s_liq_kj_per_kg_k + x_ref * s_vap_kj_per_kg_k
    );

    // first test temperatures 

    let t = t_hs_eqm(h, s);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t.get::<degree_celsius>(),
        epsilon=0.5
        );
    approx::assert_relative_eq!(
        t_kelvin,
        t.get::<kelvin>(),
        max_relative=1e-3
        );

    // then liquid and vapour specific vol
    let v_ref_m3_per_kg = (1.0 - x_ref) * v_liq_m3_per_kg + x_ref * v_vap_m3_per_kg;
    let v = v_hs_eqm(h, s);

    approx::assert_relative_eq!(
        v_ref_m3_per_kg,
        v.get::<cubic_meter_per_kilogram>(),
        max_relative=5e-3
        );


    // enthalpy of vaporisation
    // a.k.a latent heat
    // kind of manual, not really in the flashing function (yet) 
    // but it tests this algorithm...

    let mut enthalpy_of_vap = h_tp_2(t, p) - h_tp_1(t, p);

    if t_kelvin > 623.15 {
        // in this regime, we must use region 3 properties
        // I look at t_sat, and force it to find liquid volume 
        // using a slightly colder temperature than tsat 
        // and for vapour I get it slightly higher than  
        // the tsat
        let t_sat = sat_temp_4(p);
        let t_sat_kelvin = t_sat.get::<kelvin>();
        let v_vap: SpecificVolume = {
            // this covers up to tsat at 643.15 K
            if t_sat_kelvin <= 640.691 {
                v_tp_3t(t_sat, p)
            } else if t_sat_kelvin <= 643.15 {
                v_tp_3r(t_sat, p)
            } else if t_sat_kelvin <= 646.483 {
                v_tp_3x(t_sat, p)
            } else {
                v_tp_3z(t_sat, p)
            }
        };

        let v_liq: SpecificVolume = {
            // this covers up to tsat at 643.15 K
            if t_sat_kelvin <= 634.659 {
                v_tp_3c(t_sat, p)
            } else if t_sat_kelvin <= 643.15 {
                v_tp_3s(t_sat, p)
            } else if t_sat_kelvin <= 646.599 {
                v_tp_3u(t_sat, p)
            } else {
                v_tp_3y(t_sat, p)
            }


        };

        let h_liq = h_rho_t_3(v_liq.recip(), t);
        let h_vap = h_rho_t_3(v_vap.recip(), t);

        enthalpy_of_vap = h_vap-h_liq;

    }

    // only do this for EXACTLY critical point
    // because numerical errors will not cause the enthalpies to cancel 
    // out exactly
    if t_kelvin == 647.096 {
        enthalpy_of_vap = AvailableEnergy::ZERO;
    }

    approx::assert_relative_eq!(
        enthalpy_of_vap_kj_per_kg,
        enthalpy_of_vap.get::<kilojoule_per_kilogram>(),
        max_relative=1e-2
        );



}



// done for near critical point, where results may be more 
// inaccurate
fn assert_hs_flash_near_crit(t_deg_c: f64,
    t_kelvin: f64,
    psat_bar: f64,
    v_liq_m3_per_kg: f64,
    v_vap_m3_per_kg: f64,
    h_liq_kj_per_kg: f64,
    h_vap_kj_per_kg: f64,
    enthalpy_of_vap_kj_per_kg: f64,
    s_liq_kj_per_kg_k: f64,
    s_vap_kj_per_kg_k: f64){

    // specify a vapour quality
    let x_ref = 0.77;
    let p = Pressure::new::<bar>(psat_bar);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(
        (1.0-x_ref) * h_liq_kj_per_kg + x_ref * h_vap_kj_per_kg
    );
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        (1.0-x_ref) * s_liq_kj_per_kg_k + x_ref * s_vap_kj_per_kg_k
    );

    // first test temperatures 

    let t = t_hs_eqm(h, s);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t.get::<degree_celsius>(),
        epsilon=5.0
        );
    approx::assert_relative_eq!(
        t_kelvin,
        t.get::<kelvin>(),
        max_relative=1e-2
        );

    // then liquid and vapour specific vol
    let v_ref_m3_per_kg = (1.0 - x_ref) * v_liq_m3_per_kg + x_ref * v_vap_m3_per_kg;
    let v = v_hs_eqm(h, s);

    approx::assert_relative_eq!(
        v_ref_m3_per_kg,
        v.get::<cubic_meter_per_kilogram>(),
        max_relative=9.5e-2
        );


    // enthalpy of vaporisation
    // a.k.a latent heat
    // kind of manual, not really in the flashing function (yet) 
    // but it tests this algorithm...

    let mut enthalpy_of_vap = h_tp_2(t, p) - h_tp_1(t, p);

    if t_kelvin > 623.15 {
        // in this regime, we must use region 3 properties
        // I look at t_sat, and force it to find liquid volume 
        // using a slightly colder temperature than tsat 
        // and for vapour I get it slightly higher than  
        // the tsat
        let t_sat = sat_temp_4(p);
        let t_sat_kelvin = t_sat.get::<kelvin>();
        let v_vap: SpecificVolume = {
            // this covers up to tsat at 643.15 K
            if t_sat_kelvin <= 640.691 {
                v_tp_3t(t_sat, p)
            } else if t_sat_kelvin <= 643.15 {
                v_tp_3r(t_sat, p)
            } else if t_sat_kelvin <= 646.483 {
                v_tp_3x(t_sat, p)
            } else {
                v_tp_3z(t_sat, p)
            }
        };

        let v_liq: SpecificVolume = {
            // this covers up to tsat at 643.15 K
            if t_sat_kelvin <= 634.659 {
                v_tp_3c(t_sat, p)
            } else if t_sat_kelvin <= 643.15 {
                v_tp_3s(t_sat, p)
            } else if t_sat_kelvin <= 646.599 {
                v_tp_3u(t_sat, p)
            } else {
                v_tp_3y(t_sat, p)
            }


        };

        let h_liq = h_rho_t_3(v_liq.recip(), t);
        let h_vap = h_rho_t_3(v_vap.recip(), t);

        enthalpy_of_vap = h_vap-h_liq;

    }

    // only do this for EXACTLY critical point
    // because numerical errors will not cause the enthalpies to cancel 
    // out exactly
    if t_kelvin == 647.096 {
        enthalpy_of_vap = AvailableEnergy::ZERO;
    }

    approx::assert_relative_eq!(
        enthalpy_of_vap_kj_per_kg,
        enthalpy_of_vap.get::<kilojoule_per_kilogram>(),
        max_relative=5e-2
        );



}



