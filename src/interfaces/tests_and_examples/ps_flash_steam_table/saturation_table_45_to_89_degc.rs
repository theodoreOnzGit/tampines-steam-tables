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
pub fn saturation_table_45_to_89_degc(){

    //[t_deg_c,t_kelvin,psat_bar,v_liq_m3_per_kg,v_vap_m3_per_kg,h_liq_kj_per_kg,h_vap_kj_per_kg,enthalpy_of_vap,s_liq_kj_per_kg_k,s_vap_kj_per_kg_k],
    let steam_table: Vec<[f64; 10]> =
        vec![
        [45.0,318.15,0.0959439,0.00100991,15.2534,188.437,2582.45,2394.02,0.63862,8.1634],
        [46.0,319.15,0.100988,0.00101034,14.5355,192.617,2584.23,2391.61,0.65174,8.1454],
        [47.0,320.15,0.106259,0.00101078,13.8562,196.796,2586.0,2389.21,0.66481,8.1276],
        [48.0,321.15,0.111764,0.00101123,13.2132,200.976,2587.77,2386.8,0.67785,8.1099],
        [49.0,322.15,0.117512,0.00101168,12.6045,205.156,2589.54,2384.39,0.69084,8.0923],
        [50.0,323.15,0.123513,0.00101214,12.0279,209.336,2591.31,2381.97,0.70379,8.0749],
        [51.0,324.15,0.129774,0.0010126,11.4815,213.517,2593.08,2379.56,0.71671,8.0576],
        [52.0,325.15,0.136305,0.00101308,10.9637,217.697,2594.84,2377.14,0.72958,8.0405],
        [53.0,326.15,0.143116,0.00101356,10.4726,221.878,2596.6,2374.72,0.74242,8.0235],
        [54.0,327.15,0.150215,0.00101404,10.0069,226.059,2598.35,2372.3,0.75522,8.0066],
        [55.0,328.15,0.157614,0.00101454,9.56492,230.241,2600.11,2369.87,0.76798,7.9899],
        [56.0,329.15,0.165322,0.00101504,9.14543,234.423,2601.86,2367.44,0.7807,7.9733],
        [57.0,330.15,0.17335,0.00101555,8.74712,238.605,2603.61,2365.01,0.79339,7.9568],
        [58.0,331.15,0.181708,0.00101606,8.36879,242.788,2605.36,2362.57,0.80603,7.9405],
        [59.0,332.15,0.190407,0.00101658,8.00932,246.971,2607.1,2360.13,0.81864,7.9243],
        [60.0,333.15,0.199458,0.00101711,7.66766,251.154,2608.85,2357.69,0.83122,7.9082],
        [61.0,334.15,0.208873,0.00101765,7.34281,255.338,2610.58,2355.25,0.84375,7.8922],
        [62.0,335.15,0.218664,0.00101819,7.03384,259.523,2612.32,2352.8,0.85625,7.8764],
        [63.0,336.15,0.228842,0.00101874,6.7399,263.708,2614.05,2350.35,0.86872,7.8607],
        [64.0,337.15,0.239421,0.00101929,6.46015,267.893,2615.78,2347.89,0.88115,7.8451],
        [65.0,338.15,0.250411,0.00101985,6.19383,272.079,2617.51,2345.43,0.89354,7.8296],
        [66.0,339.15,0.261827,0.00102042,5.94021,276.266,2619.23,2342.97,0.9059,7.8142],
        [67.0,340.15,0.27368,0.001021,5.69861,280.453,2620.96,2340.5,0.91823,7.799],
        [68.0,341.15,0.285986,0.00102158,5.4684,284.641,2622.67,2338.03,0.93052,7.7839],
        [69.0,342.15,0.298756,0.00102216,5.24896,288.829,2624.39,2335.56,0.94277,7.7689],
        [70.0,343.15,0.312006,0.00102276,5.03973,293.018,2626.1,2333.08,0.95499,7.754],
        [71.0,344.15,0.32575,0.00102336,4.84018,297.208,2627.81,2330.6,0.96718,7.7392],
        [72.0,345.15,0.340001,0.00102396,4.6498,301.398,2629.51,2328.11,0.97933,7.7245],
        [73.0,346.15,0.354775,0.00102458,4.46812,305.589,2631.21,2325.62,0.99146,7.71],
        [74.0,347.15,0.370088,0.0010252,4.29469,309.781,2632.91,2323.13,1.0035,7.6955],
        [75.0,348.15,0.385954,0.00102582,4.12908,313.974,2634.6,2320.63,1.0156,7.6812],
        [76.0,349.15,0.402389,0.00102645,3.9709,318.167,2636.29,2318.13,1.0276,7.6669],
        [77.0,350.15,0.419409,0.00102709,3.81978,322.361,2637.98,2315.62,1.0396,7.6528],
        [78.0,351.15,0.437031,0.00102773,3.67535,326.556,2639.66,2313.11,1.0516,7.6388],
        [79.0,352.15,0.455271,0.00102838,3.53729,330.752,2641.34,2310.59,1.0635,7.6248],
        [80.0,353.15,0.474147,0.00102904,3.40527,334.949,2643.01,2308.07,1.0754,7.611],
        [81.0,354.15,0.493676,0.0010297,3.27899,339.146,2644.68,2305.54,1.0873,7.5973],
        [82.0,355.15,0.513875,0.00103037,3.15818,343.345,2646.35,2303.01,1.0991,7.5837],
        [83.0,356.15,0.534762,0.00103105,3.04257,347.544,2648.01,2300.47,1.1109,7.5701],
        [84.0,357.15,0.556355,0.00103173,2.9319,351.745,2649.67,2297.93,1.1227,7.5567],
        [85.0,358.15,0.578675,0.00103242,2.82593,355.946,2651.33,2295.38,1.1344,7.5434],
        [86.0,359.15,0.601738,0.00103311,2.72445,360.148,2652.98,2292.83,1.1461,7.5301],
        [87.0,360.15,0.625565,0.00103381,2.62722,364.352,2654.62,2290.27,1.1578,7.517],
        [88.0,361.15,0.650174,0.00103451,2.53406,368.556,2656.26,2287.7,1.1694,7.5039],
        [89.0,362.15,0.675587,0.00103522,2.44476,372.762,2657.9,2285.14,1.1811,7.4909],
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
        max_relative=1e-4
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
