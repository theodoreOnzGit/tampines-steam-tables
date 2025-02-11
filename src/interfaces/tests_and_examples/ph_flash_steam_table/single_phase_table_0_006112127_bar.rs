use uom::si::{dynamic_viscosity::micropascal_second, thermal_conductivity::milliwatt_per_meter_kelvin};
use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::{dynamic_viscosity::mu_ph_eqm, interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, lambda_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm}};

/// single phase table (see page 192)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
#[test]
pub fn single_phase_table_0_to_240_degc_except_triple_pt(){

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
            assert_ph_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}
#[test]
pub fn single_phase_table_250_to_800_degc(){

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

    // assert temp first to within 0.025 mk 
    let temp_tol_millikelvin = 25.0;

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
        max_relative=1e-4
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
        max_relative=1e-4
        );

    // dynamic_viscosity
    //
    let eta_micropascal_second_test = mu_ph_eqm(p, h)
        .get::<micropascal_second>();
    approx::assert_relative_eq!(
        eta_micropascal_second,
        eta_micropascal_second_test,
        max_relative=2e-2
        );

    // thermal thermal conductivity
    let lambda_test_milliwatt_per_meter_kelvin = 
        lambda_ph_eqm(p, h).get::<milliwatt_per_meter_kelvin>();
    approx::assert_relative_eq!(
        lambda_milliwatt_per_meter_kelvin,
        lambda_test_milliwatt_per_meter_kelvin,
        max_relative=1e-2
        );

}



