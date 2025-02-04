use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::{ph_flash_eqm::x_ph_flash, pt_flash_eqm::{cp_tp_eqm_two_phase, h_tp_eqm_two_phase, kappa_tp_eqm_two_phase, s_tp_eqm_two_phase, v_tp_eqm_two_phase, w_tp_eqm_two_phase}};

/// single phase table (see page 192)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
#[test]
pub fn single_phase_table_0_to_240_degc_except_triple_pt_0_00611_bar(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[0.006112127,0.0,0.00100021,-0.04159,-0.0001545,4.2199,1402.3,3216538.0,1792.0,555.57],
        //[0.006112127,0.0,206.14,2500.89,9.1558,1.8882,408.88,1.3269,8.9455,16.76],
        [0.006112127,2.0,207.657,2504.66,9.1695,1.8822,410.5,1.3277,9.0033,16.889],
        [0.006112127,4.0,209.173,2508.42,9.1831,1.878,412.08,1.3282,9.0617,17.019],
        [0.006112127,6.0,210.688,2512.18,9.1966,1.875,413.63,1.3286,9.1207,17.15],
        [0.006112127,8.0,212.203,2515.92,9.21,1.873,415.15,1.3288,9.1803,17.282],
        [0.006112127,10.0,213.717,2519.67,9.2233,1.8716,416.65,1.3289,9.2404,17.414],
        [0.006112127,12.0,215.231,2523.41,9.2364,1.8706,418.13,1.329,9.3011,17.548],
        [0.006112127,14.0,216.744,2527.15,9.2495,1.87,419.6,1.329,9.3622,17.682],
        [0.006112127,16.0,218.258,2530.89,9.2625,1.8696,421.06,1.329,9.4239,17.817],
        [0.006112127,18.0,219.771,2534.63,9.2754,1.8694,422.51,1.329,9.4861,17.953],
        [0.006112127,20.0,221.284,2538.37,9.2882,1.8693,423.95,1.3289,9.5488,18.089],
        [0.006112127,25.0,225.065,2547.71,9.3198,1.8695,427.53,1.3287,9.7075,18.434],
        [0.006112127,30.0,228.846,2557.06,9.3509,1.8701,431.07,1.3285,9.8689,18.784],
        [0.006112127,35.0,232.626,2566.42,9.3815,1.8708,434.57,1.3282,10.033,19.138],
        [0.006112127,40.0,236.406,2575.77,9.4116,1.8718,438.03,1.3279,10.199,19.498],
        [0.006112127,45.0,240.185,2585.13,9.4413,1.8728,441.47,1.3276,10.368,19.862],
        [0.006112127,50.0,243.964,2594.5,9.4705,1.874,444.87,1.3272,10.538,20.231],
        [0.006112127,55.0,247.743,2603.87,9.4993,1.8753,448.25,1.3269,10.711,20.604],
        [0.006112127,60.0,251.521,2613.25,9.5276,1.8767,451.59,1.3265,10.885,20.982],
        [0.006112127,65.0,255.299,2622.64,9.5556,1.8781,454.9,1.3262,11.061,21.364],
        [0.006112127,70.0,259.077,2632.03,9.5832,1.8797,458.19,1.3258,11.239,21.751],
        [0.006112127,75.0,262.855,2641.44,9.6104,1.8813,461.44,1.3254,11.419,22.142],
        [0.006112127,80.0,266.632,2650.85,9.6372,1.8831,464.67,1.3249,11.599,22.537],
        [0.006112127,85.0,270.409,2660.27,9.6637,1.8849,467.88,1.3245,11.782,22.936],
        [0.006112127,90.0,274.186,2669.7,9.6898,1.8867,471.05,1.324,11.965,23.339],
        [0.006112127,95.0,277.963,2679.14,9.7157,1.8887,474.2,1.3236,12.15,23.747],
        [0.006112127,100.0,281.74,2688.58,9.7411,1.8907,477.33,1.3231,12.336,24.158],
        [0.006112127,110.0,289.294,2707.51,9.7912,1.8949,483.51,1.3221,12.712,24.993],
        [0.006112127,120.0,296.847,2726.48,9.8401,1.8993,489.59,1.3211,13.092,25.843],
        [0.006112127,130.0,304.399,2745.5,9.8878,1.9039,495.58,1.3201,13.475,26.708],
        [0.006112127,140.0,311.952,2764.56,9.9346,1.9087,501.49,1.319,13.862,27.587],
        [0.006112127,150.0,319.504,2783.67,9.9803,1.9137,507.31,1.3179,14.252,28.481],
        [0.006112127,160.0,327.056,2802.84,10.025,1.9188,513.05,1.3168,14.645,29.388],
        [0.006112127,170.0,334.608,2822.05,10.069,1.924,518.72,1.3156,15.04,30.309],
        [0.006112127,180.0,342.16,2841.32,10.112,1.9294,524.31,1.3145,15.438,31.242],
        [0.006112127,190.0,349.712,2860.64,10.154,1.9348,529.83,1.3133,15.838,32.189],
        [0.006112127,200.0,357.264,2880.01,10.195,1.9404,535.28,1.3121,16.24,33.148],
        [0.006112127,210.0,364.816,2899.45,10.236,1.946,540.66,1.3109,16.643,34.119],
        [0.006112127,220.0,372.367,2918.93,10.276,1.9517,545.98,1.3097,17.048,35.102],
        [0.006112127,230.0,379.919,2938.48,10.315,1.9575,551.23,1.3085,17.454,36.096],
        [0.006112127,240.0,387.47,2958.08,10.354,1.9633,556.43,1.3073,17.862,37.102],
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
            assert_pt_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}
#[test]
pub fn single_phase_table_250_to_800_degc_0_00611_bar(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        [0.006112127,250.0,395.022,2977.75,10.392,1.9692,561.56,1.3061,18.27,38.119],
        [0.006112127,260.0,402.573,2997.47,10.429,1.9752,566.64,1.3049,18.68,39.147],
        [0.006112127,270.0,410.125,3017.25,10.466,1.9812,571.67,1.3037,19.09,40.185],
        [0.006112127,280.0,417.676,3037.09,10.502,1.9873,576.64,1.3025,19.501,41.234],
        [0.006112127,290.0,425.227,3057.0,10.538,1.9934,581.56,1.3013,19.913,42.292],
        [0.006112127,300.0,432.779,3076.96,10.573,1.9996,586.43,1.3001,20.325,43.361],
        [0.006112127,310.0,440.33,3096.99,10.608,2.0058,591.25,1.2989,20.737,44.439],
        [0.006112127,320.0,447.881,3117.08,10.642,2.012,596.02,1.2977,21.15,45.527],
        [0.006112127,330.0,455.432,3137.23,10.675,2.0183,600.74,1.2965,21.563,46.623],
        [0.006112127,340.0,462.984,3157.44,10.709,2.0247,605.42,1.2953,21.977,47.729],
        [0.006112127,350.0,470.535,3177.72,10.741,2.031,610.06,1.2941,22.39,48.844],
        [0.006112127,360.0,478.086,3198.07,10.774,2.0374,614.65,1.2929,22.803,49.967],
        [0.006112127,370.0,485.637,3218.47,10.806,2.0439,619.2,1.2917,23.217,51.099],
        [0.006112127,380.0,493.188,3238.94,10.837,2.0504,623.7,1.2905,23.63,52.239],
        [0.006112127,390.0,500.74,3259.48,10.869,2.0569,628.17,1.2893,24.043,53.388],
        [0.006112127,400.0,508.291,3280.08,10.899,2.0635,632.6,1.2881,24.456,54.544],
        [0.006112127,410.0,515.842,3300.75,10.93,2.0701,636.99,1.2869,24.868,55.708],
        [0.006112127,420.0,523.393,3321.48,10.96,2.0767,641.34,1.2858,25.281,56.88],
        [0.006112127,430.0,530.944,3342.28,10.99,2.0834,645.65,1.2846,25.692,58.059],
        [0.006112127,440.0,538.495,3363.15,11.019,2.0901,649.93,1.2834,26.104,59.245],
        [0.006112127,450.0,546.046,3384.08,11.048,2.0968,654.18,1.2822,26.515,60.438],
        [0.006112127,460.0,553.598,3405.09,11.077,2.1036,658.39,1.2811,26.925,61.639],
        [0.006112127,470.0,561.149,3426.16,11.106,2.1103,662.56,1.2799,27.335,62.846],
        [0.006112127,480.0,568.7,3447.29,11.134,2.1172,666.7,1.2788,27.745,64.06],
        [0.006112127,490.0,576.251,3468.5,11.162,2.124,670.81,1.2776,28.154,65.281],
        [0.006112127,500.0,583.802,3489.77,11.19,2.1309,674.89,1.2765,28.562,66.508],
        [0.006112127,510.0,591.353,3511.12,11.217,2.1378,678.94,1.2753,28.969,67.741],
        [0.006112127,520.0,598.904,3532.53,11.244,2.1447,682.96,1.2742,29.376,68.981],
        [0.006112127,530.0,606.455,3554.01,11.271,2.1517,686.95,1.2731,29.783,70.226],
        [0.006112127,540.0,614.006,3575.56,11.298,2.1586,690.9,1.272,30.188,71.478],
        [0.006112127,550.0,621.557,3597.18,11.324,2.1656,694.83,1.2708,30.593,72.735],
        [0.006112127,560.0,629.108,3618.88,11.351,2.1726,698.74,1.2697,30.997,73.998],
        [0.006112127,570.0,636.66,3640.64,11.376,2.1796,702.61,1.2686,31.4,75.266],
        [0.006112127,580.0,644.211,3662.47,11.402,2.1867,706.46,1.2675,31.802,76.54],
        [0.006112127,590.0,651.762,3684.37,11.428,2.1937,710.29,1.2664,32.204,77.819],
        [0.006112127,600.0,659.313,3706.34,11.453,2.2008,714.08,1.2654,32.605,79.104],
        [0.006112127,650.0,697.068,3817.27,11.577,2.2362,732.71,1.2601,34.595,85.6],
        [0.006112127,700.0,734.823,3929.96,11.695,2.2716,750.77,1.255,36.564,92.212],
        [0.006112127,750.0,772.578,4044.43,11.81,2.307,768.31,1.2501,38.508,98.926],
        [0.006112127,800.0,810.333,4160.66,11.921,2.3423,785.38,1.2454,40.428,105.73],

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
            assert_pt_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}

/// single phase table (see page 192)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_0_to_240_degc_except_triple_pt_0_1_bar(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[0.1,0.0,0.0010002,-0.03202,-0.0001539,4.2199,1402.3,196604.0,1792.0,555.58],
        [0.1,2.0,0.0010001,8.40098,0.030607,4.2133,1412.1,199379.0,1673.7,560.6],
        [0.1,4.0,0.00100007,16.8219,0.061101,4.2078,1421.5,202050.0,1567.4,565.4],
        [0.1,6.0,0.0010001,25.2326,0.091339,4.2031,1430.5,204616.0,1471.6,570.01],
        [0.1,8.0,0.00100019,33.6348,0.12133,4.1991,1439.2,207075.0,1384.8,574.45],
        [0.1,10.0,0.00100034,42.0296,0.15108,4.1958,1447.4,209428.0,1306.0,578.72],
        [0.1,12.0,0.00100054,50.4183,0.18061,4.1929,1455.3,211674.0,1234.1,582.83],
        [0.1,14.0,0.0010008,58.8017,0.2099,4.1905,1462.8,213814.0,1168.4,586.81],
        [0.1,16.0,0.0010011,67.1805,0.23898,4.1884,1470.0,215848.0,1108.1,590.65],
        [0.1,18.0,0.00100145,75.5555,0.26785,4.1866,1476.8,217779.0,1052.7,594.36],
        [0.1,20.0,0.00100184,83.9271,0.2965,4.1851,1483.3,219606.0,1001.6,597.96],
        [0.1,25.0,0.001003,104.845,0.36725,4.1822,1498.0,223737.0,890.04,606.46],
        [0.1,30.0,0.00100441,125.75,0.43679,4.1803,1510.8,227261.0,797.22,614.35],
        [0.1,35.0,0.00100604,146.649,0.50517,4.1792,1521.8,230209.0,719.12,621.66],
        [0.1,40.0,0.00100788,167.543,0.57243,4.1788,1531.2,232610.0,652.72,628.45],
        [0.1,45.0,0.00100991,188.438,0.63862,4.179,1538.9,234495.0,595.76,634.75],
        [0.1,45.8075,0.00101026,191.812,0.64922,4.1791,1540.0,234753.0,587.32,635.72],
        [0.1,45.8075,14.6706,2583.89,8.1489,1.9413,440.51,1.3227,10.377,19.942],
        [0.1,50.0,14.8674,2591.99,8.1741,1.9272,443.67,1.324,10.521,20.255],
        [0.1,55.0,15.1015,2601.6,8.2037,1.9178,447.26,1.3246,10.694,20.631],
        [0.1,60.0,15.3353,2611.18,8.2326,1.9124,450.74,1.3248,10.87,21.011],
        [0.1,65.0,15.5687,2620.73,8.2611,1.9091,454.14,1.3247,11.047,21.396],
        [0.1,70.0,15.8018,2630.27,8.2891,1.9071,457.5,1.3246,11.226,21.784],
        [0.1,75.0,16.0347,2639.8,8.3167,1.9058,460.81,1.3243,11.406,22.176],
        [0.1,80.0,16.2674,2649.33,8.3438,1.9051,464.09,1.324,11.587,22.572],
        [0.1,85.0,16.4999,2658.86,8.3706,1.9048,467.33,1.3236,11.77,22.972],
        [0.1,90.0,16.7323,2668.38,8.397,1.9048,470.55,1.3233,11.955,23.376],
        [0.1,95.0,16.9646,2677.9,8.4231,1.9051,473.73,1.3229,12.14,23.784],
        [0.1,100.0,17.1967,2687.43,8.4488,1.9057,476.89,1.3225,12.327,24.196],
        [0.1,110.0,17.6607,2706.5,8.4992,1.9074,483.12,1.3216,12.703,25.03],
        [0.1,120.0,18.1243,2725.58,8.5484,1.9098,489.25,1.3207,13.084,25.88],
        [0.1,130.0,18.5876,2744.69,8.5964,1.9128,495.28,1.3197,13.468,26.744],
        [0.1,140.0,19.0507,2763.84,8.6433,1.9163,501.22,1.3187,13.856,27.622],
        [0.1,150.0,19.5136,2783.02,8.6892,1.9201,507.07,1.3176,14.247,28.515],
        [0.1,160.0,19.9763,2802.24,8.734,1.9243,512.83,1.3166,14.64,29.42],
        [0.1,170.0,20.4389,2821.51,8.778,1.9288,518.52,1.3154,15.036,30.34],
        [0.1,180.0,20.9013,2840.82,8.8211,1.9336,524.13,1.3143,15.434,31.272],
        [0.1,190.0,21.3637,2860.18,8.8634,1.9385,529.66,1.3132,15.834,32.217],
        [0.1,200.0,21.826,2879.59,8.9048,1.9436,535.13,1.312,16.236,33.175],
        [0.1,210.0,22.2882,2899.05,8.9455,1.9489,540.52,1.3108,16.64,34.145],
        [0.1,220.0,22.7504,2918.57,8.9855,1.9543,545.85,1.3097,17.045,35.126],
        [0.1,230.0,23.2124,2938.14,9.0248,1.9598,551.12,1.3085,17.452,36.119],
        [0.1,240.0,23.6745,2957.76,9.0634,1.9654,556.32,1.3073,17.859,37.124],
        ];

        for dataset in steam_table {
            let p_bar = dataset[0] as f64;
            let t_deg_c = dataset[1] as f64;
            let v_m3_per_kg = dataset[2] as f64;
            let h_kj_per_kg = dataset[3] as f64;
            let s_kj_per_kg_k = dataset[4] as f64;
            let cp_kj_per_kg_k = dataset[5] as f64;
            let w_m_per_s = dataset[6] as f64;
            let kappa_dimensionless = dataset[7] as f64;
            let eta_micropascal_second = dataset[8] as f64;
            let lambda_milliwatt_per_meter_kelvin = dataset[9] as f64;
            assert_pt_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}
#[test]
pub fn single_phase_table_250_to_800_degc_0_1_bar(){

    let steam_table: Vec<[f64; 10]> =
        vec![

        [0.1,250.0,24.1365,2977.45,9.1014,1.9711,561.46,1.3061,18.268,38.14],
        [0.1,260.0,24.5985,2997.19,9.1388,1.9769,566.55,1.3049,18.678,39.17],
        [0.1,270.0,25.0604,3016.98,9.1756,1.9827,571.58,1.3037,19.089,40.2],
        [0.1,280.0,25.5223,3036.84,9.2118,1.9887,576.56,1.3025,19.5,41.25],
        [0.1,290.0,25.9842,3056.76,9.2475,1.9947,581.48,1.3013,19.912,42.31],
        [0.1,300.0,26.446,3076.73,9.2827,2.0007,586.36,1.3001,20.324,43.38],
        [0.1,310.0,26.9078,3096.77,9.3173,2.0069,591.18,1.2989,20.736,44.45],
        [0.1,320.0,27.3696,3116.87,9.3515,2.013,595.96,1.2977,21.149,45.54],
        [0.1,330.0,27.8314,3137.03,9.3852,2.0192,600.68,1.2965,21.563,46.64],
        [0.1,340.0,28.2932,3157.26,9.4185,2.0255,605.37,1.2953,21.976,47.74],
        [0.1,350.0,28.755,3177.54,9.4513,2.0318,610.0,1.2941,22.389,48.86],
        [0.1,360.0,29.2167,3197.89,9.4837,2.0382,614.6,1.2929,22.803,49.98],
        [0.1,370.0,29.6785,3218.31,9.5157,2.0446,619.15,1.2917,23.216,51.11],
        [0.1,380.0,30.1402,3238.79,9.5473,2.051,623.66,1.2905,23.63,52.25],
        [0.1,390.0,30.6019,3259.33,9.5785,2.0575,628.13,1.2893,24.043,53.4],
        [0.1,400.0,31.0636,3279.94,9.6093,2.0641,632.56,1.2881,24.456,54.55],
        [0.1,410.0,31.5253,3300.61,9.6398,2.0706,636.95,1.2869,24.868,55.72],
        [0.1,420.0,31.987,3321.35,9.6699,2.0772,641.3,1.2857,25.28,56.89],
        [0.1,430.0,32.4486,3342.15,9.6997,2.0839,645.62,1.2846,25.692,58.07],
        [0.1,440.0,32.9103,3363.03,9.7292,2.0905,649.9,1.2834,26.104,59.25],
        [0.1,450.0,33.372,3383.96,9.7584,2.0972,654.15,1.2822,26.515,60.45],
        [0.1,460.0,33.8336,3404.97,9.7872,2.104,658.36,1.2811,26.925,61.65],
        [0.1,470.0,34.2953,3426.04,9.8158,2.1107,662.53,1.2799,27.335,62.85],
        [0.1,480.0,34.7569,3447.19,9.844,2.1175,666.68,1.2788,27.745,64.07],
        [0.1,490.0,35.2185,3468.4,9.872,2.1244,670.79,1.2776,28.154,65.29],
        [0.1,500.0,35.6802,3489.67,9.8997,2.1312,674.87,1.2765,28.562,66.52],
        [0.1,510.0,36.1418,3511.02,9.9271,2.1381,678.92,1.2753,28.97,67.75],
        [0.1,520.0,36.6034,3532.44,9.9543,2.145,682.94,1.2742,29.377,68.99],
        [0.1,530.0,37.065,3553.92,9.9812,2.152,686.92,1.2731,29.783,70.23],
        [0.1,540.0,37.5267,3575.48,10.008,2.1589,690.88,1.272,30.188,71.48],
        [0.1,550.0,37.9883,3597.1,10.034,2.1659,694.82,1.2708,30.593,72.74],
        [0.1,560.0,38.4499,3618.79,10.061,2.1729,698.72,1.2697,30.997,74.0],
        [0.1,570.0,38.9115,3640.56,10.086,2.1799,702.6,1.2686,31.4,75.27],
        [0.1,580.0,39.3731,3662.39,10.112,2.1869,706.45,1.2675,31.803,76.55],
        [0.1,590.0,39.8347,3684.3,10.138,2.194,710.27,1.2664,32.204,77.83],
        [0.1,600.0,40.2963,3706.27,10.163,2.201,714.07,1.2654,32.605,79.11],
        [0.1,650.0,42.6042,3817.2,10.287,2.2364,732.7,1.2601,34.596,85.61],
        [0.1,700.0,44.9121,3929.91,10.405,2.2718,750.76,1.255,36.564,92.22],
        [0.1,750.0,47.2199,4044.38,10.52,2.3071,768.31,1.2501,38.508,98.93],
        [0.1,800.0,49.5278,4160.62,10.631,2.3424,785.38,1.2454,40.429,105.74],


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
            assert_pt_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}

/// single phase table (see page 192)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_2_to_800_degc_except_triple_pt_1_bar(){

    let steam_table: Vec<[f64; 10]> =
        vec![

        //[1.0,0.0,0.00100016,0.05966,-0.0001478,4.2194,1402.4,19665.0,1791.8,555.65],
        [1.0,2.0,0.00100006,8.49179,0.03061,4.2129,1412.2,19943.0,1673.5,560.66],
        [1.0,4.0,0.00100003,16.9119,0.061101,4.2074,1421.6,20210.0,1567.3,565.47],
        [1.0,6.0,0.00100006,25.3219,0.091336,4.2027,1430.7,20467.0,1471.5,570.08],
        [1.0,8.0,0.00100015,33.7233,0.12133,4.1988,1439.3,20713.0,1384.7,574.51],
        [1.0,10.0,0.0010003,42.1174,0.15108,4.1955,1447.6,20948.0,1305.9,578.78],
        [1.0,12.0,0.0010005,50.5054,0.1806,4.1926,1455.4,21172.0,1234.0,582.89],
        [1.0,14.0,0.00100076,58.8881,0.20989,4.1902,1463.0,21386.0,1168.3,586.86],
        [1.0,16.0,0.00100106,67.2664,0.23897,4.1881,1470.1,21590.0,1108.1,590.7],
        [1.0,18.0,0.00100141,75.6407,0.26783,4.1863,1476.9,21783.0,1052.7,594.42],
        [1.0,20.0,0.0010018,84.0118,0.29648,4.1848,1483.4,21966.0,1001.6,598.01],
        [1.0,25.0,0.00100296,104.928,0.36723,4.1819,1498.2,22379.0,890.02,606.52],
        [1.0,30.0,0.00100437,125.833,0.43676,4.18,1511.0,22731.0,797.22,614.39],
        [1.0,35.0,0.001006,146.73,0.50513,4.179,1522.0,23026.0,719.13,621.71],
        [1.0,40.0,0.00100784,167.623,0.57239,4.1786,1531.3,23266.0,652.73,628.49],
        [1.0,45.0,0.00100987,188.516,0.63859,4.1788,1539.0,23455.0,595.77,634.8],
        [1.0,50.0,0.0010121,209.412,0.70375,4.1796,1545.3,23595.0,546.52,640.64],
        [1.0,55.0,0.0010145,230.313,0.76794,4.1809,1550.2,23689.0,503.63,646.04],
        [1.0,60.0,0.00101708,251.222,0.83117,4.1828,1553.9,23739.0,466.04,651.02],
        [1.0,65.0,0.00101982,272.141,0.8935,4.1852,1556.3,23749.0,432.91,655.59],
        [1.0,70.0,0.00102273,293.074,0.95495,4.1881,1557.6,23722.0,403.56,659.78],
        [1.0,75.0,0.00102579,314.023,1.0156,4.1915,1557.8,23658.0,377.42,663.58],
        [1.0,80.0,0.00102902,334.991,1.0754,4.1955,1557.1,23561.0,354.06,667.01],
        [1.0,85.0,0.00103239,355.979,1.1344,4.2,1555.4,23432.0,333.08,670.08],
        [1.0,90.0,0.00103593,376.992,1.1926,4.205,1552.8,23275.0,314.18,672.8],
        [1.0,95.0,0.00103962,398.03,1.2502,4.2106,1549.3,23090.0,297.09,675.17],
        [1.0,99.6,0.00104315,417.436,1.3026,4.2161,1545.5,22896.0,282.75,677.07],
        [1.0,99.6,1.69402,2674.95,7.3588,2.0759,472.05,1.3154,12.218,24.532],
        [1.0,100.0,1.69596,2675.77,7.361,2.0741,472.34,1.3155,12.234,24.564],
        [1.0,110.0,1.74482,2696.32,7.4154,2.0399,479.27,1.3165,12.62,25.397],
        [1.0,120.0,1.79324,2716.61,7.4676,2.0187,485.89,1.3166,13.009,26.24],
        [1.0,130.0,1.84132,2736.72,7.5181,2.0039,492.31,1.3163,13.401,27.096],
        [1.0,140.0,1.88913,2756.7,7.5671,1.9933,498.57,1.3158,13.796,27.963],
        [1.0,150.0,1.93673,2776.59,7.6147,1.9857,504.7,1.3152,14.192,28.843],
        [1.0,160.0,1.98414,2796.42,7.661,1.9805,510.7,1.3145,14.591,29.736],
        [1.0,170.0,2.0314,2816.21,7.7062,1.9772,516.59,1.3137,14.992,30.642],
        [1.0,180.0,2.07853,2835.97,7.7503,1.9755,522.38,1.3129,15.394,31.56],
        [1.0,190.0,2.12556,2855.72,7.7934,1.9751,528.07,1.3119,15.798,32.491],
        [1.0,200.0,2.17249,2875.48,7.8356,1.9757,533.67,1.311,16.204,33.436],
        [1.0,210.0,2.21935,2895.24,7.8769,1.9772,539.19,1.3099,16.611,34.392],
        [1.0,220.0,2.26614,2915.02,7.9174,1.9793,544.62,1.3089,17.019,35.361],
        [1.0,230.0,2.31287,2934.83,7.9572,1.9821,549.98,1.3078,17.428,36.342],
        [1.0,240.0,2.35955,2954.66,7.9962,1.9854,555.27,1.3067,17.838,37.335],
        [1.0,250.0,2.40619,2974.54,8.0346,1.9891,560.49,1.3056,18.249,38.34],
        [1.0,260.0,2.45279,2994.45,8.0723,1.9932,565.65,1.3045,18.661,39.356],
        [1.0,270.0,2.49935,3014.4,8.1094,1.9975,570.74,1.3033,19.073,40.383],
        [1.0,280.0,2.54588,3034.4,8.1458,2.0022,575.77,1.3022,19.486,41.421],
        [1.0,290.0,2.59239,3054.45,8.1818,2.007,580.75,1.301,19.899,42.47],
        [1.0,300.0,2.63887,3074.54,8.2171,2.0121,585.67,1.2998,20.313,43.53],
        [1.0,310.0,2.68533,3094.69,8.252,2.0173,590.54,1.2987,20.727,44.599],
        [1.0,320.0,2.73176,3114.89,8.2863,2.0227,595.35,1.2975,21.141,45.679],
        [1.0,330.0,2.77818,3135.14,8.3202,2.0282,600.11,1.2963,21.555,46.768],
        [1.0,340.0,2.82458,3155.45,8.3536,2.0338,604.83,1.2951,21.969,47.867],
        [1.0,350.0,2.87097,3175.82,8.3865,2.0396,609.5,1.2939,22.384,48.975],
        [1.0,360.0,2.91735,3196.24,8.419,2.0454,614.12,1.2928,22.798,50.093],
        [1.0,370.0,2.96371,3216.73,8.4511,2.0514,618.7,1.2916,23.212,51.219],
        [1.0,380.0,3.01006,3237.27,8.4828,2.0574,623.23,1.2904,23.626,52.354],
        [1.0,390.0,3.05639,3257.87,8.5141,2.0635,627.73,1.2892,24.04,53.498],
        [1.0,400.0,3.10272,3278.54,8.5451,2.0697,632.18,1.2881,24.453,54.649],
        [1.0,410.0,3.14904,3299.27,8.5756,2.0759,636.59,1.2869,24.866,55.809],
        [1.0,420.0,3.19535,3320.06,8.6059,2.0822,640.96,1.2857,25.279,56.977],
        [1.0,430.0,3.24165,3340.91,8.6357,2.0886,645.3,1.2845,25.692,58.153],
        [1.0,440.0,3.28795,3361.83,8.6653,2.095,649.59,1.2834,26.103,59.336],
        [1.0,450.0,3.33424,3382.81,8.6945,2.1015,653.85,1.2822,26.515,60.527],
        [1.0,460.0,3.38052,3403.86,8.7234,2.108,658.08,1.2811,26.926,61.725],
        [1.0,470.0,3.42679,3424.97,8.752,2.1146,662.27,1.2799,27.336,62.93],
        [1.0,480.0,3.47306,3446.15,8.7803,2.1212,666.43,1.2788,27.746,64.142],
        [1.0,490.0,3.51932,3467.4,8.8083,2.1279,670.55,1.2776,28.155,65.361],
        [1.0,500.0,3.56558,3488.71,8.8361,2.1345,674.64,1.2765,28.564,66.586],
        [1.0,510.0,3.61184,3510.09,8.8635,2.1413,678.7,1.2753,28.971,67.818],
        [1.0,520.0,3.65809,3531.53,8.8907,2.148,682.73,1.2742,29.379,69.056],
        [1.0,530.0,3.70433,3553.05,8.9177,2.1548,686.73,1.2731,29.785,70.3],
        [1.0,540.0,3.75057,3574.63,8.9444,2.1617,690.7,1.272,30.191,71.551],
        [1.0,550.0,3.79681,3596.28,8.9709,2.1685,694.64,1.2709,30.596,72.807],
        [1.0,560.0,3.84304,3618.0,8.9971,2.1754,698.55,1.2698,31.0,74.069],
        [1.0,570.0,3.88928,3639.79,9.0231,2.1823,702.44,1.2687,31.403,75.337],
        [1.0,580.0,3.9355,3661.65,9.0489,2.1892,706.29,1.2676,31.806,76.611],
        [1.0,590.0,3.98173,3683.58,9.0744,2.1962,710.12,1.2665,32.207,77.889],
        [1.0,600.0,4.02795,3705.57,9.0998,2.2031,713.93,1.2654,32.608,79.173],
        [1.0,650.0,4.25902,3816.6,9.2234,2.2381,732.59,1.2601,34.6,85.67],
        [1.0,700.0,4.49004,3929.38,9.3424,2.2732,750.68,1.255,36.568,92.282],
        [1.0,750.0,4.72101,4043.92,9.4571,2.3083,768.24,1.2502,38.512,98.999],
        [1.0,800.0,4.95196,4160.21,9.5681,2.3434,785.34,1.2455,40.433,105.81],

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
            assert_pt_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }
}


fn assert_pt_flash(
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
    let t = ThermodynamicTemperature::new::<degree_celsius>(t_deg_c);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(h_kj_per_kg);
    let x = x_ph_flash(p, h_ref);

    dbg!(&x);

    // assert h
    let h_test = h_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        h_kj_per_kg,
        h_test.get::<kilojoule_per_kilogram>(),
        max_relative=1e-4
        );

    // assert volume to within 0.01%  (that's the tolerable error for 
    // backward eqn)
    let v_test = v_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        v_m3_per_kg,
        v_test.get::<cubic_meter_per_kilogram>(),
        max_relative=1e-4
        );

    // now entropy 
    let s_test = s_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        s_kj_per_kg_k,
        s_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-3
        );

    // cp 
    let cp_test = cp_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        cp_kj_per_kg_k,
        cp_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-3
        );
    // w 
    let w_test = w_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=1e-4
        );

    // kappa
    let kappa_test = kappa_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=1e-3
        );

    // eta and lambda tbd



}



