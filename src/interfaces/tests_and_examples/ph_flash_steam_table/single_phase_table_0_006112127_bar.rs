use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm};
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;

#[test]
pub fn reminder_to_do_viscosity_and_thermal_conductivity(){
    todo!()
}
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

    // eta and lambda tbd



}



