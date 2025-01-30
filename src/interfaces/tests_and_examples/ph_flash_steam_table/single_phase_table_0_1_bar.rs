use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm};

/// single phase table (see page 192)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
#[test]
pub fn single_phase_table_0_to_240_degc_except_triple_pt(){

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
            assert_ph_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}
pub fn single_phase_table_250_to_800_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![

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
        max_relative=5e-4
        );

    // eta and lambda tbd



}



