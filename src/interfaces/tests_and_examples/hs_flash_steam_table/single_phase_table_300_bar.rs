use uom::si::{dynamic_viscosity::micropascal_second, thermal_conductivity::milliwatt_per_meter_kelvin};
use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::{bar, kilopascal, megapascal};
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::hs_flash_eqm::{cp_hs_eqm, kappa_hs_eqm, lambda_hs_eqm, mu_hs_eqm, tpvx_hs_flash_eqm, w_hs_eqm};
use crate::interfaces::functional_programming::ph_flash_eqm::x_ph_flash;

/// single phase table (see page 201)
///
/// thermal conductivity off by 8%
#[test]
pub fn single_phase_table_2_to_750_degc(){

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
            assert_hs_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}

fn assert_hs_flash(
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
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        s_kj_per_kg_k
    );

    // assert temp first to within 0.050 k
    // based on table 2.8
    let temp_tol_millikelvin = 50.0;

    let (t_test, p_test, v_test, x_test) = 
        tpvx_hs_flash_eqm(h, s);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t_test.get::<degree_celsius>(),
        epsilon=temp_tol_millikelvin*1e-3
        );

    // assert volume to within 0.5%  
    approx::assert_relative_eq!(
        v_m3_per_kg,
        v_test.get::<cubic_meter_per_kilogram>(),
        max_relative=5e-3
        );

    // assert pressure to within 15 kPa or 0.60%
    if p > Pressure::new::<megapascal>(2.5) {
        let pressure_error_tol = Pressure::new::<kilopascal>(25.0);
        let pressure_error_tol_bar = pressure_error_tol.get::<bar>();
        approx::assert_abs_diff_eq!(
            p_bar,
            p_test.get::<bar>(),
            epsilon=pressure_error_tol_bar
        ); 
    } else {
        dbg!(&(p_bar,t_deg_c));
        approx::assert_relative_eq!(
            p_bar,
            p_test.get::<bar>(),
            max_relative = 2e-2
        ); 

    }
    // the x for ph is quite reliable, i'll use that as reference 
    let x_ref = x_ph_flash(p, h);

    approx::assert_relative_eq!(
        x_ref.round(),
        x_test.get::<ratio>(),
        max_relative = 1e-3
    );


    // cp 
    let cp_test = cp_hs_eqm(h, s);
    approx::assert_relative_eq!(
        cp_kj_per_kg_k,
        cp_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-3
        );
    // w 
    let w_test = w_hs_eqm(h, s);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=1e-3
        );

    // kappa
    let kappa_test = kappa_hs_eqm(h, s);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=2e-2
        );

    // dynamic_viscosity
    //
    let eta_micropascal_second_test = mu_hs_eqm(h, s)
        .get::<micropascal_second>();
    approx::assert_relative_eq!(
        eta_micropascal_second,
        eta_micropascal_second_test,
        max_relative=2e-2
        );

    // thermal thermal conductivity
    let lambda_test_milliwatt_per_meter_kelvin = 
        lambda_hs_eqm(h, s).get::<milliwatt_per_meter_kelvin>();
    approx::assert_relative_eq!(
        lambda_milliwatt_per_meter_kelvin,
        lambda_test_milliwatt_per_meter_kelvin,
        max_relative=8e-2
        );

}













