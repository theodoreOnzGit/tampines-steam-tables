//! Zaloudek validation: multiphase IN-DOME stagnation critical flow.
//! Validates `get_critical_pressure_and_mass_flux_ph_vle_dome` on stagnation
//! states that lie inside the p-h VLE dome (filtered via region check).

use uom::si::f64::*;
use uom::si::area::square_foot;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::mass_rate::pound_per_second;
use uom::si::pressure::{kilopascal, pound_force_per_square_inch};

use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
use crate::interfaces::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::prelude::functional_programming::pt_flash_eqm::FwdEqnRegion;
use crate::steam_turbine_equations::choked_flow::get_stagnation_conditions_from_throat_ph;
use crate::steam_turbine_equations::choked_flow::get_critical_pressure_and_mass_flux_ph_vle_dome;


/// Validates the in-dome critical-flow solver
/// `get_critical_pressure_and_mass_flux_ph_vle_dome` against Zaloudek data.
///
/// Method:
///   1. Build the throat state from (p_crit_ref, x_t).
///   2. Backward map (throat -> stagnation) to recover the stagnation (p0, h0).
///   3. Filter: if the stagnation (p0, h0) lies OUTSIDE the VLE dome
///      (i.e. not Region4), skip it — it belongs to a different bucket
///      (subcooled liquid or superheated/supercritical vapour) and is not
///      what this solver is for.
///   4. For in-dome stagnation points, run the in-dome solver and assert the
///      recovered critical pressure matches Zaloudek's p_crit and the mass
///      flux matches Zaloudek's G_crit.
fn validate_zaloudek_curve_in_dome(
    x_t: f64,
    data: &[(f64, f64, f64)],
    critical_pressure_tolerance: f64,
    mass_flux_log_tolerance: f64,
) {
    let ref_vol = TampinesSteamTableCV::get_ref_vol();

    for &(p_psia, g_expected_val, _h0_expected_val) in data {
        let p_throat_critical_ref = Pressure::new::<pound_force_per_square_inch>(p_psia);
        let g_expected = MassRate::new::<pound_per_second>(g_expected_val)
            / Area::new::<square_foot>(1.0);

        // throat state from saturation pressure + quality
        let state_t = TampinesSteamTableCV::new_from_sat_pressure_quality(
            p_throat_critical_ref, x_t, ref_vol);
        let h_t = state_t.get_specific_enthalpy();

        // backward map: throat -> stagnation
        let (p0, h0, _g_throat) =
            get_stagnation_conditions_from_throat_ph(p_throat_critical_ref, h_t);

        // filter: only test stagnation points that lie inside the dome
        if ph_flash_region(p0, h0) != FwdEqnRegion::Region4 {
            eprintln!(
                "skip p={p_psia} psia, x_t={x_t}: stagnation lies outside dome \
                 ({:?})",
                ph_flash_region(p0, h0)
            );
            continue;
        }

        // in-dome solver under test
        let (p_crit_calc, g_calc) =
            get_critical_pressure_and_mass_flux_ph_vle_dome(p0, h0);

        dbg!(&(p_psia, x_t,
               p_crit_calc.get::<kilopascal>(),
               p_throat_critical_ref.get::<kilopascal>(),
               g_calc.get::<kilogram_per_square_meter_second>(),
               g_expected.get::<kilogram_per_square_meter_second>()));

        approx::assert_relative_eq!(
            p_crit_calc.get::<kilopascal>(),
            p_throat_critical_ref.get::<kilopascal>(),
            max_relative = critical_pressure_tolerance,
        );

        approx::assert_relative_eq!(
            g_calc.get::<kilogram_per_square_meter_second>().log10(),
            g_expected.get::<kilogram_per_square_meter_second>().log10(),
            max_relative = mass_flux_log_tolerance,
        );
    }
}


#[test]
fn quality_0_50_in_dome(){
    // throat quality x_t = 0.50 (50%)
    let data: Vec<(f64, f64, f64)> = vec![
        (5.0,    26.3991,   650.2463),
        (10.0,   50.4252,   667.9803),
        (15.0,   73.7254,   679.8030),
        (20.0,   96.3178,   697.5369),
        (30.0,   144.8425,  715.2709),
        (50.0,   247.2153,  735.9606),
        (75.0,   356.3976,  753.6946),
        (100.0,  465.6123,  765.5172),
        (150.0,  700.1867,  789.1626),
        (200.0,  940.8566,  800.9852),
        (300.0,  1337.4357, 821.6749),
        (500.0,  2157.8051, 848.2759),
        (750.0,  3290.8767, 871.9212),
        (1000.0, 4484.6769, 883.7438),
        (1500.0, 6198.1302, 901.4778),
        (2000.0, 8446.5673, 913.3005),
        (3000.0, 12006.8680,913.3005),
    ];
    validate_zaloudek_curve_in_dome(0.50, &data, 0.005, 0.05);
}
