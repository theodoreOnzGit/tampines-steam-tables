use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm};

/// single phase table (see page 201)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_2_to_800_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[6.0,0.0,0.000999902,0.56879,-0.0001144,4.2169,1403.2,3282.0,1790.6,556.03],
        [6.0,2.0,0.000999806,8.9961,0.030625,4.2105,1413.0,3328.3,1672.6,561.03],
        [6.0,4.0,0.000999778,17.4117,0.0611,4.2052,1422.4,3372.9,1566.5,565.82],
        [6.0,6.0,0.000999814,25.8173,0.09132,4.2006,1431.4,3415.7,1470.8,570.42],
        [6.0,8.0,0.000999909,34.2147,0.12129,4.1968,1440.1,3456.7,1384.2,574.84],
        [6.0,10.0,0.00100006,42.605,0.15103,4.1936,1448.3,3495.9,1305.5,579.1],
        [6.0,12.0,0.00100027,50.9892,0.18054,4.1908,1456.2,3533.4,1233.7,583.21],
        [6.0,14.0,0.00100052,59.3684,0.20982,4.1884,1463.8,3569.1,1168.0,587.17],
        [6.0,16.0,0.00100082,67.7432,0.23889,4.1864,1470.9,3603.0,1107.8,591.01],
        [6.0,18.0,0.00100117,76.1143,0.26774,4.1847,1477.7,3635.2,1052.5,594.71],
        [6.0,20.0,0.00100157,84.4823,0.29638,4.1832,1484.2,3665.7,1001.4,598.3],
        [6.0,25.0,0.00100274,105.391,0.3671,4.1805,1499.0,3734.6,889.95,606.8],
        [6.0,30.0,0.00100414,126.288,0.43661,4.1787,1511.8,3793.4,797.21,614.67],
        [6.0,35.0,0.00100578,147.179,0.50496,4.1777,1522.8,3842.6,719.16,621.97],
        [6.0,40.0,0.00100762,168.066,0.5722,4.1773,1532.1,3882.7,652.79,628.76],
        [6.0,45.0,0.00100965,188.953,0.63837,4.1776,1539.9,3914.2,595.86,635.06],
        [6.0,50.0,0.00101188,209.843,0.70352,4.1784,1546.2,3937.6,546.62,640.9],
        [6.0,55.0,0.00101428,230.738,0.76769,4.1798,1551.1,3953.3,503.74,646.3],
        [6.0,60.0,0.00101685,251.642,0.83091,4.1817,1554.7,3961.8,466.16,651.28],
        [6.0,65.0,0.00101959,272.556,0.89322,4.1841,1557.2,3963.6,433.04,655.85],
        [6.0,70.0,0.0010225,293.483,0.95465,4.187,1558.5,3959.0,403.69,660.04],
        [6.0,75.0,0.00102556,314.426,1.0152,4.1905,1558.7,3948.5,377.56,663.84],
        [6.0,80.0,0.00102878,335.388,1.075,4.1944,1558.0,3932.4,354.19,667.28],
        [6.0,85.0,0.00103215,356.372,1.134,4.1989,1556.3,3911.1,333.22,670.35],
        [6.0,90.0,0.00103568,377.378,1.1923,4.2039,1553.7,3884.9,314.32,673.07],
        [6.0,95.0,0.00103937,398.412,1.2498,4.2094,1550.3,3854.2,297.22,675.45],
        [6.0,100.0,0.0010432,419.474,1.3066,4.2155,1546.1,3819.1,281.72,677.5],
        [6.0,110.0,0.00105134,461.696,1.4183,4.2293,1535.4,3737.3,254.73,680.62],
        [6.0,120.0,0.0010601,504.067,1.5275,4.2453,1521.9,3641.5,232.14,682.48],
        [6.0,130.0,0.00106952,546.611,1.6343,4.2639,1505.8,3533.5,213.03,683.16],
        [6.0,140.0,0.00107961,589.355,1.739,4.2853,1487.2,3414.6,196.7,682.68],
        [6.0,150.0,0.00109042,632.328,1.8418,4.3099,1466.3,3286.3,182.64,681.1],
        [6.0,158.832,0.00110061,670.501,1.9311,4.3345,1445.9,3166.0,171.77,679.02],
        [6.0,158.832,0.315575,2756.14,6.7592,2.48,495.88,1.2987,14.264,31.552],
        [6.0,160.0,0.316666,2759.02,6.7658,2.4597,496.87,1.2994,14.314,31.645],
        [6.0,170.0,0.325824,2782.97,6.8205,2.3427,504.69,1.3029,14.743,32.45],
        [6.0,180.0,0.334744,2806.04,6.872,2.2759,511.81,1.3042,15.171,33.273],
        [6.0,190.0,0.343496,2828.55,6.9211,2.2295,518.55,1.3047,15.597,34.111],
        [6.0,200.0,0.352116,2850.66,6.9684,2.1942,525.04,1.3048,16.023,34.967],
        [6.0,210.0,0.360628,2872.46,7.0139,2.1666,531.32,1.3047,16.448,35.839],
        [6.0,220.0,0.36905,2894.01,7.0581,2.1446,537.43,1.3044,16.873,36.727],
        [6.0,230.0,0.377395,2915.37,7.101,2.1273,543.38,1.3039,17.296,37.632],
        [6.0,240.0,0.385675,2936.57,7.1427,2.1136,549.19,1.3034,17.72,38.553],
        [6.0,250.0,0.393899,2957.65,7.1834,2.103,554.88,1.3027,18.143,39.49],
        [6.0,260.0,0.402075,2978.64,7.2231,2.0949,560.45,1.302,18.566,40.442],
        [6.0,270.0,0.410209,2999.56,7.262,2.0888,565.91,1.3012,18.988,41.41],
        [6.0,280.0,0.418307,3020.42,7.3001,2.0845,571.28,1.3003,19.41,42.392],
        [6.0,290.0,0.426372,3041.25,7.3374,2.0817,576.56,1.2994,19.831,43.389],
        [6.0,300.0,0.434409,3062.06,7.374,2.08,581.76,1.2985,20.253,44.399],
        [6.0,310.0,0.44242,3082.86,7.41,2.0794,586.88,1.2975,20.673,45.424],
        [6.0,320.0,0.450409,3103.65,7.4453,2.0797,591.92,1.2965,21.094,46.461],
        [6.0,330.0,0.458378,3124.45,7.4801,2.0807,596.89,1.2954,21.514,47.511],
        [6.0,340.0,0.466328,3145.27,7.5143,2.0823,601.8,1.2944,21.934,48.573],
        [6.0,350.0,0.474262,3166.1,7.548,2.0845,606.65,1.2933,22.353,49.648],
        [6.0,360.0,0.482181,3186.96,7.5812,2.0872,611.43,1.2922,22.772,50.734],
        [6.0,370.0,0.490086,3207.85,7.614,2.0903,616.16,1.2911,23.19,51.831],
        [6.0,380.0,0.497979,3228.77,7.6463,2.0938,620.84,1.29,23.608,52.939],
        [6.0,390.0,0.50586,3249.72,7.6781,2.0976,625.46,1.2889,24.025,54.058],
        [6.0,400.0,0.513731,3270.72,7.7095,2.1017,630.04,1.2878,24.442,55.187],
        [6.0,410.0,0.521592,3291.76,7.7405,2.106,634.56,1.2867,24.858,56.326],
        [6.0,420.0,0.529444,3312.84,7.7712,2.1106,639.04,1.2855,25.273,57.475],
        [6.0,430.0,0.537287,3333.97,7.8014,2.1154,643.47,1.2844,25.688,58.633],
        [6.0,440.0,0.545123,3355.15,7.8314,2.1203,647.86,1.2833,26.102,59.801],
        [6.0,450.0,0.552952,3376.38,7.8609,2.1255,652.21,1.2822,26.516,60.977],
        [6.0,460.0,0.560774,3397.66,7.8901,2.1307,656.52,1.281,26.929,62.162],
        [6.0,470.0,0.568589,3418.99,7.919,2.1362,660.79,1.2799,27.341,63.355],
        [6.0,480.0,0.576399,3440.38,7.9476,2.1417,665.02,1.2788,27.752,64.556],
        [6.0,490.0,0.584203,3461.83,7.9759,2.1473,669.21,1.2776,28.163,65.765],
        [6.0,500.0,0.592002,3483.33,8.0039,2.1531,673.37,1.2765,28.573,66.982],
        [6.0,510.0,0.599796,3504.89,8.0316,2.159,677.49,1.2754,28.982,68.206],
        [6.0,520.0,0.607586,3526.51,8.0591,2.1649,681.58,1.2743,29.39,69.437],
        [6.0,530.0,0.615371,3548.19,8.0862,2.1709,685.63,1.2732,29.798,70.676],
        [6.0,540.0,0.623153,3569.93,8.1131,2.177,689.66,1.2721,30.205,71.921],
        [6.0,550.0,0.63093,3591.73,8.1398,2.1832,693.65,1.271,30.61,73.173],
        [6.0,560.0,0.638704,3613.59,8.1662,2.1895,697.61,1.2699,31.015,74.431],
        [6.0,570.0,0.646475,3635.52,8.1923,2.1958,701.54,1.2688,31.42,75.696],
        [6.0,580.0,0.654242,3657.51,8.2183,2.2021,705.44,1.2677,31.823,76.967],
        [6.0,590.0,0.662006,3679.56,8.2439,2.2085,709.31,1.2667,32.225,78.243],
        [6.0,600.0,0.669768,3701.68,8.2694,2.215,713.16,1.2656,32.627,79.526],
        [6.0,650.0,0.708536,3813.24,8.3937,2.2477,731.99,1.2604,34.62,86.021],
        [6.0,700.0,0.747253,3926.46,8.5131,2.2811,750.22,1.2553,36.59,92.64],
        [6.0,750.0,0.785928,4041.36,8.6282,2.3148,767.91,1.2505,38.535,99.369],
        [6.0,800.0,0.824571,4157.95,8.7395,2.3488,785.1,1.2459,40.455,106.19],
        ];

        for dataset in steam_table {
            let p_bar = dataset[0];
            let t_deg_c = dataset[1];
            let v_m3_per_kg = dataset[2];
            let h_kj_per_kg = dataset[3];
            let s_kj_per_kg_k = dataset[4];
            let cp_kj_per_kg_k = dataset[5];
            let w_m_per_s = dataset[6];
            let kappa_dimensionless = dataset[7];
            let eta_micropascal_second = dataset[8];
            let lambda_milliwatt_per_meter_kelvin = dataset[9];
            assert_ph_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }
}

fn assert_ph_flash(
    p_bar: f64,
    t_deg_c: f64,
    v_m3_per_kg: f64,
    h_kj_per_kg: f64,
    s_kj_per_kg_k: f64,
    cp_kj_per_kg_k: f64,
    w_m_per_s: f64,
    kappa_dimensionless: f64,
    eta_micropascal_second: f64,
    lambda_milliwatt_per_meter_kelvin: f64,
){
    let p = Pressure::new::<bar>(p_bar);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(h_kj_per_kg);

    // assert temp first to within 0.035 mk 
    let temp_tol_millikelvin = 35.0;

    let t_test = t_ph_eqm(p, h);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t_test.get::<degree_celsius>(),
        epsilon=temp_tol_millikelvin*1e-3
        );

    // assert volume to within 0.01%  (that's the tolerable error for 
    // backward eqn)
    let v_test = v_ph_eqm(p, h);
    approx::assert_relative_eq!(
        v_m3_per_kg,
        v_test.get::<cubic_meter_per_kilogram>(),
        max_relative=1e-4
        );

    // now entropy 
    let s_test = s_ph_eqm(p, h);
    approx::assert_relative_eq!(
        s_kj_per_kg_k,
        s_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=7e-3
        );

    // cp 
    let cp_test = cp_ph_eqm(p, h);
    approx::assert_relative_eq!(
        cp_kj_per_kg_k,
        cp_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-4
        );
    // w 
    let w_test = w_ph_eqm(p, h);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=1e-4
        );

    // kappa
    let kappa_test = kappa_ph_eqm(p, h);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=6e-3
        );

    // eta and lambda tbd



}



