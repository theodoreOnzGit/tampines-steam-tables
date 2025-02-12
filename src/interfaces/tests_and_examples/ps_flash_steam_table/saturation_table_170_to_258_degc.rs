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
pub fn saturation_table_170_to_258_degc(){

    //[t_deg_c,t_kelvin,psat_bar,v_liq_m3_per_kg,v_vap_m3_per_kg,h_liq_kj_per_kg,h_vap_kj_per_kg,enthalpy_of_vap,s_liq_kj_per_kg_k,s_vap_kj_per_kg_k],
    let steam_table: Vec<[f64; 10]> =
        vec![
        [170.0,443.15,7.92053,0.00111426,0.242616,719.206,2767.89,2048.69,2.0419,6.6649],
        [172.0,445.15,8.31077,0.00111682,0.231785,727.973,2769.85,2041.88,2.0616,6.6485],
        [174.0,447.15,8.71606,0.00111941,0.221528,736.755,2771.77,2035.01,2.0811,6.6322],
        [176.0,449.15,9.13681,0.00112203,0.21181,745.551,2773.63,2028.08,2.1007,6.6161],
        [178.0,451.15,9.57343,0.00112469,0.202598,754.362,2775.45,2021.09,2.1201,6.6],
        [180.0,453.15,10.0263,0.00112739,0.193862,763.188,2777.22,2014.03,2.1395,6.5841],
        [182.0,455.15,10.496,0.00113012,0.185572,772.03,2778.94,2006.91,2.1589,6.5682],
        [184.0,457.15,10.9827,0.00113289,0.177703,780.889,2780.61,1999.72,2.1782,6.5525],
        [186.0,459.15,11.4871,0.0011357,0.170229,789.764,2782.23,1992.47,2.1974,6.5369],
        [188.0,461.15,12.0094,0.00113855,0.163127,798.656,2783.8,1985.14,2.2166,6.5214],
        [190.0,463.15,12.5502,0.00114144,0.156377,807.566,2785.31,1977.74,2.2358,6.506],
        [192.0,465.15,13.1099,0.00114437,0.149957,816.494,2786.77,1970.28,2.2549,6.4907],
        [194.0,467.15,13.6889,0.00114734,0.143848,825.44,2788.18,1962.74,2.2739,6.4755],
        [196.0,469.15,14.2877,0.00115036,0.138034,834.405,2789.53,1955.12,2.2929,6.4603],
        [198.0,471.15,14.9069,0.00115341,0.132497,843.389,2790.82,1947.44,2.3119,6.4453],
        [200.0,473.15,15.5467,0.00115651,0.127222,852.393,2792.06,1939.67,2.3308,6.4303],
        [202.0,475.15,16.2078,0.00115966,0.122195,861.417,2793.24,1931.82,2.3497,6.4154],
        [204.0,477.15,16.8906,0.00116285,0.117402,870.463,2794.36,1923.9,2.3685,6.4006],
        [206.0,479.15,17.5955,0.00116609,0.11283,879.529,2795.42,1915.89,2.3873,6.3858],
        [208.0,481.15,18.3231,0.00116937,0.108467,888.618,2796.42,1907.8,2.406,6.3711],
        [210.0,483.15,19.0739,0.00117271,0.104302,897.729,2797.35,1899.62,2.4248,6.3565],
        [212.0,485.15,19.8483,0.00117609,0.100325,906.863,2798.22,1891.36,2.4434,6.342],
        [214.0,487.15,20.647,0.00117953,0.0965249,916.021,2799.03,1883.01,2.4621,6.3275],
        [216.0,489.15,21.4702,0.00118302,0.0928934,925.203,2799.77,1874.57,2.4807,6.313],
        [218.0,491.15,22.3187,0.00118656,0.0894214,934.409,2800.45,1866.04,2.4993,6.2986],
        [220.0,493.15,23.1929,0.00119016,0.0861007,943.642,2801.05,1857.41,2.5178,6.2842],
        [222.0,495.15,24.0933,0.00119381,0.0829236,952.9,2801.59,1848.69,2.5363,6.2699],
        [224.0,497.15,25.0205,0.00119752,0.0798826,962.185,2802.05,1839.87,2.5548,6.2557],
        [226.0,499.15,25.9749,0.00120129,0.076971,971.498,2802.45,1830.95,2.5733,6.2414],
        [228.0,501.15,26.9572,0.00120512,0.0741823,980.839,2802.76,1821.93,2.5917,6.2272],
        [230.0,503.15,27.9679,0.00120901,0.0715102,990.21,2803.01,1812.8,2.6102,6.2131],
        [232.0,505.15,29.0075,0.00121297,0.0689492,999.609,2803.18,1803.57,2.6285,6.1989],
        [234.0,507.15,30.0767,0.00121699,0.0664936,1009.04,2803.27,1794.23,2.6469,6.1848],
        [236.0,509.15,31.1758,0.00122108,0.0641385,1018.5,2803.28,1784.78,2.6653,6.1707],
        [238.0,511.15,32.3056,0.00122523,0.0618788,1028.0,2803.21,1775.22,2.6836,6.1566],
        [240.0,513.15,33.4665,0.00122946,0.0597101,1037.52,2803.06,1765.54,2.7019,6.1425],
        [242.0,515.15,34.6592,0.00123376,0.057628,1047.08,2802.82,1755.74,2.7203,6.1285],
        [244.0,517.15,35.8843,0.00123814,0.0556284,1056.68,2802.5,1745.82,2.7385,6.1144],
        [246.0,519.15,37.1423,0.00124259,0.0537073,1066.31,2802.1,1735.78,2.7568,6.1003],
        [248.0,521.15,38.4338,0.00124712,0.0518612,1075.98,2801.6,1725.62,2.7751,6.0863],
        [250.0,523.15,39.7594,0.00125174,0.0500866,1085.69,2801.01,1715.33,2.7934,6.0722],
        [252.0,525.15,41.1197,0.00125644,0.0483801,1095.43,2800.33,1704.9,2.8117,6.0582],
        [254.0,527.15,42.5154,0.00126122,0.0467386,1105.22,2799.56,1694.34,2.8299,6.0441],
        [256.0,529.15,43.9471,0.0012661,0.0451592,1115.04,2798.69,1683.64,2.8482,6.03],
        [258.0,531.15,45.4153,0.00127107,0.043639,1124.91,2797.71,1672.8,2.8664,6.0158],
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
