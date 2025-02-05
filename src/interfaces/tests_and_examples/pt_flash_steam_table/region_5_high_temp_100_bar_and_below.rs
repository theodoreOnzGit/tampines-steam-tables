use uom::si::{available_energy::kilojoule_per_kilogram, temperature_interval};
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;

use crate::{interfaces::functional_programming::{ph_flash_eqm::x_ph_flash, pt_flash_eqm::{cp_tp_eqm_two_phase, h_tp_eqm_two_phase, kappa_tp_eqm_two_phase, s_tp_eqm_two_phase, v_tp_eqm_two_phase}}, region_4_vap_liq_equilibrium::sat_temp_4};

#[test]
pub fn single_phase_table_0_00611_bar(){

    let steam_table: Vec<[f64; 7]> =
        vec![
        [800.0,0.006112127,810.333,4160.66,11.921,2.3423,1.2454],
        [820.0,0.006112127,825.435,4207.69,11.964,2.3573,1.2435],
        [840.0,0.006112127,840.538,4254.97,12.007,2.3707,1.2417],
        [860.0,0.006112127,855.64,4302.52,12.05,2.3841,1.2401],
        [880.0,0.006112127,870.742,4350.33,12.091,2.3975,1.2384],
        [900.0,0.006112127,885.844,4398.42,12.133,2.4109,1.2368],
        [920.0,0.006112127,900.946,4446.77,12.174,2.4243,1.2351],
        [940.0,0.006112127,916.048,4495.39,12.214,2.4376,1.2336],
        [960.0,0.006112127,931.15,4544.27,12.254,2.4508,1.232],
        [980.0,0.006112127,946.252,4593.42,12.294,2.4639,1.2305],
        [1000.0,0.006112127,961.354,4642.82,12.333,2.4769,1.229],
        [1020.0,0.006112127,976.456,4692.49,12.371,2.4898,1.2275],
        [1040.0,0.006112127,991.558,4742.42,12.41,2.5025,1.2261],
        [1060.0,0.006112127,1006.66,4792.59,12.448,2.5152,1.2247],
        [1080.0,0.006112127,1021.76,4843.02,12.485,2.5276,1.2234],
        [1100.0,0.006112127,1036.86,4893.7,12.522,2.54,1.2221],
        [1120.0,0.006112127,1051.97,4944.62,12.559,2.5521,1.2208],
        [1140.0,0.006112127,1067.07,4995.78,12.596,2.5641,1.2195],
        [1160.0,0.006112127,1082.17,5047.18,12.632,2.5759,1.2183],
        [1180.0,0.006112127,1097.27,5098.82,12.668,2.5876,1.2171],
        [1200.0,0.006112127,1112.37,5150.68,12.703,2.5991,1.2159],
        [1220.0,0.006112127,1127.48,5202.78,12.738,2.6104,1.2148],
        [1240.0,0.006112127,1142.58,5255.1,12.773,2.6215,1.2137],
        [1260.0,0.006112127,1157.68,5307.64,12.807,2.6325,1.2126],
        [1280.0,0.006112127,1172.78,5360.4,12.842,2.6432,1.2115],
        [1300.0,0.006112127,1187.88,5413.37,12.876,2.6538,1.2105],
        [1320.0,0.006112127,1202.99,5466.55,12.909,2.6642,1.2095],
        [1340.0,0.006112127,1218.09,5519.93,12.942,2.6745,1.2086],
        [1360.0,0.006112127,1233.19,5573.52,12.975,2.6845,1.2076],
        [1380.0,0.006112127,1248.29,5627.31,13.008,2.6944,1.2067],
        [1400.0,0.006112127,1263.39,5681.3,13.041,2.7041,1.2058],
        [1420.0,0.006112127,1278.5,5735.48,13.073,2.7136,1.2049],
        [1440.0,0.006112127,1293.6,5789.84,13.105,2.723,1.2041],
        [1460.0,0.006112127,1308.7,5844.39,13.136,2.7322,1.2033],
        [1480.0,0.006112127,1323.8,5899.13,13.168,2.7412,1.2025],
        [1500.0,0.006112127,1338.9,5954.04,13.199,2.7501,1.2017],
        [1550.0,0.006112127,1376.66,6092.08,13.276,2.7715,1.1998],
        [1600.0,0.006112127,1414.41,6231.18,13.351,2.7921,1.198],
        [1650.0,0.006112127,1452.17,6371.27,13.425,2.8117,1.1964],
        [1700.0,0.006112127,1489.92,6512.33,13.497,2.8305,1.1948],
        [1750.0,0.006112127,1527.68,6654.32,13.568,2.8486,1.1933],
        [1800.0,0.006112127,1565.43,6797.19,13.638,2.8661,1.1919],
        [1850.0,0.006112127,1603.19,6940.91,13.707,2.8829,1.1906],
        [1900.0,0.006112127,1640.94,7085.47,13.774,2.8992,1.1893],
        [1950.0,0.006112127,1678.7,7230.83,13.84,2.9152,1.1881],
        [2000.0,0.006112127,1716.45,7376.98,13.905,2.9307,1.1869],

        ];

        for dataset in steam_table {
            let t_deg_c = dataset[0];
            let p_bar = dataset[1];
            let v_m3_per_kg = dataset[2];
            let h_kj_per_kg = dataset[3];
            let s_kj_per_kg_k = dataset[4];
            let cp_kj_per_kg_k = dataset[5];
            let kappa_dimensionless = dataset[6];
            assert_pt_flash_high_temp(t_deg_c, p_bar, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, kappa_dimensionless, 
                );
        }

}


fn assert_pt_flash_high_temp(
    t_deg_c: f64,
    p_bar: f64,
    v_m3_per_kg: f64,
    h_kj_per_kg: f64,
    s_kj_per_kg_k: f64,
    cp_kj_per_kg_k: f64,
    kappa_dimensionless: f64,
){
    let p = Pressure::new::<bar>(p_bar);
    let t = ThermodynamicTemperature::new::<degree_celsius>(t_deg_c);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(h_kj_per_kg);
    // check if t is less than or equal to tsat 
    
    let t_sat = sat_temp_4(p);
    let fifty_kelvin = TemperatureInterval::new::<temperature_interval::kelvin>(50.0);
    
    let x: f64;
    
    if t > t_sat + fifty_kelvin {
        x = 1.0;
    } else {
        x = x_ph_flash(p, h_ref);
    };


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

    // kappa
    let kappa_test = kappa_tp_eqm_two_phase(t, p, x);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=9e-3
        );

    // eta and lambda tbd



}




