use uom::si::{available_energy::kilojoule_per_kilogram, thermodynamic_temperature::kelvin};
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;

use crate::interfaces::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::interfaces::functional_programming::{ps_flash_eqm::t_ps_eqm, pt_flash_eqm::{h_tp_eqm_two_phase, s_tp_eqm_two_phase, v_tp_eqm_two_phase}};
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;

/// saturation table (see page 174)
#[test]
pub fn saturation_table_90_to_168_degc(){

    //[t_deg_c,t_kelvin,psat_bar,v_liq_m3_per_kg,v_vap_m3_per_kg,h_liq_kj_per_kg,h_vap_kj_per_kg,enthalpy_of_vap,s_liq_kj_per_kg_k,s_vap_kj_per_kg_k],
    let steam_table: Vec<[f64; 10]> =
        vec![
        [90.0,363.15,0.701824,0.00103594,2.35915,376.968,2659.53,2282.56,1.1927,7.4781],
        [91.0,364.15,0.728904,0.00103667,2.27705,381.176,2661.16,2279.98,1.2042,7.4653],
        [92.0,365.15,0.756849,0.0010374,2.1983,385.385,2662.78,2277.39,1.2158,7.4526],
        [93.0,366.15,0.785681,0.00103813,2.12275,389.595,2664.39,2274.8,1.2273,7.44],
        [94.0,367.15,0.81542,0.00103887,2.05025,393.806,2666.01,2272.2,1.2387,7.4275],
        [95.0,368.15,0.846089,0.00103962,1.98065,398.019,2667.61,2269.6,1.2502,7.415],
        [96.0,369.15,0.877711,0.00104038,1.91383,402.232,2669.22,2266.98,1.2616,7.4027],
        [97.0,370.15,0.910308,0.00104114,1.84965,406.447,2670.81,2264.37,1.273,7.3904],
        [98.0,371.15,0.943902,0.0010419,1.78801,410.663,2672.4,2261.74,1.2844,7.3782],
        [99.0,372.15,0.978518,0.00104268,1.72878,414.88,2673.99,2259.11,1.2957,7.3661],
        [100.0,373.15,1.01418,0.00104346,1.67186,419.099,2675.57,2256.47,1.307,7.3541],
        [102.0,375.15,1.08873,0.00104503,1.56454,427.541,2678.72,2251.18,1.3296,7.3303],
        [104.0,377.15,1.16776,0.00104663,1.46529,435.988,2681.84,2245.85,1.352,7.3068],
        [106.0,379.15,1.25147,0.00104826,1.37342,444.44,2684.94,2240.5,1.3743,7.2836],
        [108.0,381.15,1.34007,0.00104991,1.28831,452.899,2688.02,2235.12,1.3965,7.2607],
        [110.0,383.15,1.43376,0.00105158,1.20939,461.363,2691.07,2229.7,1.4187,7.238],
        [112.0,385.15,1.53277,0.00105328,1.13615,469.834,2694.09,2224.26,1.4407,7.2157],
        [114.0,387.15,1.63734,0.001055,1.06813,478.312,2697.09,2218.78,1.4626,7.1937],
        [116.0,389.15,1.74768,0.00105675,1.00489,486.796,2700.07,2213.27,1.4844,7.1719],
        [118.0,391.15,1.86404,0.00105853,0.94607,495.287,2703.02,2207.73,1.5062,7.1504],
        [120.0,393.15,1.98665,0.00106033,0.891304,503.785,2705.93,2202.15,1.5278,7.1291],
        [122.0,395.15,2.11578,0.00106215,0.840276,512.29,2708.82,2196.53,1.5494,7.1081],
        [124.0,397.15,2.25168,0.001064,0.792695,520.803,2711.69,2190.88,1.5708,7.0873],
        [126.0,399.15,2.3946,0.00106588,0.748294,529.323,2714.52,2185.19,1.5922,7.0668],
        [128.0,401.15,2.54481,0.00106778,0.706832,537.851,2717.32,2179.47,1.6134,7.0465],
        [130.0,403.15,2.7026,0.00106971,0.668084,546.388,2720.09,2173.7,1.6346,7.0264],
        [132.0,405.15,2.86823,0.00107167,0.631849,554.933,2722.83,2167.89,1.6557,7.0066],
        [134.0,407.15,3.04199,0.00107365,0.597939,563.486,2725.53,2162.04,1.6767,6.9869],
        [136.0,409.15,3.22417,0.00107566,0.566183,572.048,2728.2,2156.15,1.6977,6.9675],
        [138.0,411.15,3.41508,0.0010777,0.536425,580.62,2730.84,2150.22,1.7185,6.9483],
        [140.0,413.15,3.61501,0.00107976,0.508519,589.2,2733.44,2144.24,1.7393,6.9293],
        [142.0,415.15,3.82427,0.00108185,0.482334,597.79,2736.01,2138.22,1.76,6.9105],
        [144.0,417.15,4.04318,0.00108397,0.457748,606.39,2738.54,2132.15,1.7806,6.8918],
        [146.0,419.15,4.27205,0.00108612,0.434648,615.0,2741.04,2126.04,1.8011,6.8734],
        [148.0,421.15,4.51122,0.0010883,0.412931,623.621,2743.5,2119.88,1.8216,6.8551],
        [150.0,423.15,4.76101,0.0010905,0.392502,632.252,2745.92,2113.67,1.842,6.837],
        [152.0,425.15,5.02177,0.00109274,0.373273,640.893,2748.3,2107.41,1.8623,6.8191],
        [154.0,427.15,5.29383,0.00109501,0.355162,649.546,2750.64,2101.1,1.8825,6.8014],
        [156.0,429.15,5.57755,0.0010973,0.338095,658.211,2752.95,2094.74,1.9027,6.7838],
        [158.0,431.15,5.87329,0.00109963,0.322002,666.887,2755.21,2088.32,1.9228,6.7664],
        [160.0,433.15,6.18139,0.00110199,0.306818,675.575,2757.43,2081.86,1.9428,6.7491],
        [162.0,435.15,6.50224,0.00110438,0.292486,684.275,2759.61,2075.33,1.9627,6.732],
        [164.0,437.15,6.83619,0.0011068,0.278948,692.988,2761.75,2068.76,1.9826,6.715],
        [166.0,439.15,7.18364,0.00110925,0.266155,701.714,2763.84,2062.13,2.0025,6.6982],
        [168.0,441.15,7.54495,0.00111174,0.254059,710.453,2765.89,2055.44,2.0222,6.6815],
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
            assert_ps_flash(t_deg_c, t_kelvin, psat_bar, 
                v_liq_m3_per_kg, v_vap_m3_per_kg, h_liq_kj_per_kg, 
                h_vap_kj_per_kg, enthalpy_of_vap_kj_per_kg, 
                s_liq_kj_per_kg_k, s_vap_kj_per_kg_k);
        }

}



fn assert_ps_flash(t_deg_c: f64,
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
    let x_ref = 0.3;
    let p = Pressure::new::<bar>(psat_bar);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        (1.0-x_ref) * s_liq_kj_per_kg_k + x_ref * s_vap_kj_per_kg_k);

    // first test temperatures 

    let t = t_ps_eqm(p, s);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t.get::<degree_celsius>(),
        epsilon=1e-3
        );
    approx::assert_relative_eq!(
        t_kelvin,
        t.get::<kelvin>(),
        max_relative=1e-5
        );

    // then liquid and vapour specific vol
    let v_ref_m3_per_kg = (1.0 - x_ref) * v_liq_m3_per_kg + x_ref * v_vap_m3_per_kg;
    let v = v_ps_eqm(p, s);

    approx::assert_relative_eq!(
        v_ref_m3_per_kg,
        v.get::<cubic_meter_per_kilogram>(),
        max_relative=1e-4
        );


    let h_liq = h_tp_1(t, p);
    let h_vap = h_tp_2(t, p);

    approx::assert_relative_eq!(
        h_liq_kj_per_kg,
        h_liq.get::<kilojoule_per_kilogram>(),
        max_relative=1e-5
        );
    approx::assert_relative_eq!(
        h_vap_kj_per_kg,
        h_vap.get::<kilojoule_per_kilogram>(),
        max_relative=1e-4
        );
    // enthalpy of vaporisation
    // a.k.a latent heat
    // kind of manual, not really in the flashing function
    let enthalpy_of_vap = h_vap - h_liq;

    approx::assert_relative_eq!(
        enthalpy_of_vap_kj_per_kg,
        enthalpy_of_vap.get::<kilojoule_per_kilogram>(),
        max_relative=1e-5
        );



}
