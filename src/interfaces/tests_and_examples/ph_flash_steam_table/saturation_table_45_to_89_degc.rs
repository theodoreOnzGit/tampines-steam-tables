use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};

use crate::interfaces::functional_programming::ph_flash_eqm::{s_ph_eqm, t_ph_eqm, v_ph_eqm};
use crate::region_1_subcooled_liquid::h_tp_1;
use crate::region_2_vapour::h_tp_2;

/// saturation table (see page 174)
#[test]
pub fn saturation_table_0_to_44_degc(){

    //[t_deg_c,t_kelvin,psat_bar,v_liq_m3_per_kg,v_vap_m3_per_kg,h_liq_kj_per_kg,h_vap_kj_per_kg,enthalpy_of_vap,s_liq_kj_per_kg_k,s_vap_kj_per_kg_k],
    let steam_table: Vec<[f64; 10]> =
        vec![

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
            assert_ph_flash(t_deg_c, t_kelvin, psat_bar, 
                v_liq_m3_per_kg, v_vap_m3_per_kg, h_liq_kj_per_kg, 
                h_vap_kj_per_kg, enthalpy_of_vap_kj_per_kg, 
                s_liq_kj_per_kg_k, s_vap_kj_per_kg_k);
        }

}

fn assert_ph_flash(t_deg_c: f64,
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
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(
        (1.0-x_ref) * h_liq_kj_per_kg + x_ref * h_vap_kj_per_kg);

    // first test temperatures 

    let t = t_ph_eqm(p, h);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t.get::<degree_celsius>(),
        epsilon=1e-4
        );
    approx::assert_relative_eq!(
        t_kelvin,
        t.get::<kelvin>(),
        max_relative=1e-5
        );

    // then liquid and vapour specific vol
    let v_ref_m3_per_kg = (1.0 - x_ref) * v_liq_m3_per_kg + x_ref * v_vap_m3_per_kg;
    let v = v_ph_eqm(p, h);

    approx::assert_relative_eq!(
        v_ref_m3_per_kg,
        v.get::<cubic_meter_per_kilogram>(),
        max_relative=1e-5
        );

    // liquid and vapour s 
    let s_ref_kj_per_kg_k = 
        (1.0 - x_ref) * s_liq_kj_per_kg_k 
        + x_ref * s_vap_kj_per_kg_k;

    let s = s_ph_eqm(p, h);

    approx::assert_relative_eq!(
        s_ref_kj_per_kg_k,
        s.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-5
        );

    // enthalpy of vaporisation
    // a.k.a latent heat
    // kind of manual, not really in the flashing function
    let enthalpy_of_vap = h_tp_2(t, p) - h_tp_1(t, p);

    approx::assert_relative_eq!(
        enthalpy_of_vap_kj_per_kg,
        enthalpy_of_vap.get::<kilojoule_per_kilogram>(),
        max_relative=1e-5
        );



}
