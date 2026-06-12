//! Validation of the bubble-point locator `bubble_point_pressure_from_entropy`
//! against tabulated saturated-liquid steam-table data.
//!
//! By definition the saturated-liquid entropy at the saturation pressure,
//! `s_f(p_sat)`, must map back to `p_sat`. So for each tabulated pair
//! `(p_sat, s_f)` we feed `s_f` to the locator and assert it recovers `p_sat`.
//! This exercises both the bisection inversion and the agreement of the code's
//! saturation entropy with the steam table.
//!
//! Reference data: the `s_liq` column of the saturated steam tables in
//! `interfaces/tests_and_examples/hs_flash_steam_table/saturation_table_*`.

use uom::si::f64::*;
use uom::si::pressure::bar;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;

use crate::steam_turbine_equations::choked_flow::bubble_point_pressure_from_entropy;
use crate::steam_turbine_equations::choked_flow::dew_point_pressure_from_entropy;

#[test]
fn bubble_point_matches_saturation_table() {
    // (p_sat [bar], s_f [kJ/kg/K]) tabulated saturated-liquid data,
    // spanning 0.7 bar all the way up to the critical point (220.64 bar).
    let data: Vec<(f64, f64)> = vec![
        (0.701824, 1.1927),
        (1.01418,  1.3070),
        (1.98665,  1.5278),
        (3.61501,  1.7393),
        (6.18139,  1.9428),
        (7.54495,  2.0222),
        (9.57343,  2.1201),
        (12.0094,  2.2166),
        (14.9069,  2.3119),
        (18.3231,  2.4060),
        (22.3187,  2.4993),
        (26.9572,  2.5917),
        (32.3056,  2.6836),
        (38.4338,  2.7751),
        (45.4153,  2.8664),
        // high pressure, into the Region-3 cap (> 165.29 bar = 16.529 MPa)
        (48.4640,  2.9030),
        (64.1646,  3.0681),
        (83.4835,  3.2358),
        (98.6475,  3.3506),
        (120.505,  3.4997),
        (146.002,  3.6599),
        (165.292,  3.7783),
        (186.664,  3.9164),
        (210.434,  4.1142),
        (218.132,  4.2377),
        (220.640,  4.4120), // critical point (s_f = s_g = s_crit)
    ];

    for (psat_bar, s_f_val) in data {
        let s_f = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(s_f_val);
        let p_bubble = bubble_point_pressure_from_entropy(s_f);

        dbg!(&(s_f_val, p_bubble.get::<bar>(), psat_bar));

        approx::assert_relative_eq!(
            p_bubble.get::<bar>(),
            psat_bar,
            max_relative = 0.001,
        );
    }
}

#[test]
fn dew_point_matches_saturation_table() {
    // (p_sat [bar], s_g [kJ/kg/K]) tabulated saturated-vapour data,
    // spanning 0.7 bar all the way up to the critical point (220.64 bar).
    let data: Vec<(f64, f64)> = vec![
        (0.701824, 7.4781),
        (1.01418,  7.3541),
        (1.98665,  7.1291),
        (3.61501,  6.9293),
        (6.18139,  6.7491),
        (8.31077,  6.6485),
        (10.9827,  6.5525),
        (14.2877,  6.4603),
        (18.3231,  6.3711),
        (23.1929,  6.2842),
        (29.0075,  6.1989),
        (35.8843,  6.1144),
        (43.9471,  6.0300),
        // high pressure, into the Region-3 cap (> 165.29 bar = 16.529 MPa)
        (48.4640,  5.9875),
        (64.1646,  5.8578),
        (83.4835,  5.7215),
        (98.6475,  5.6243),
        (120.505,  5.4911),
        (146.002,  5.3359),
        (165.292,  5.2109),
        (186.664,  5.0527),
        (210.434,  4.7996),
        (218.132,  4.6299),
        (220.640,  4.4120), // critical point (s_f = s_g = s_crit)
    ];

    for (psat_bar, s_g_val) in data {
        let s_g = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(s_g_val);
        let p_dew = dew_point_pressure_from_entropy(s_g);

        dbg!(&(s_g_val, p_dew.get::<bar>(), psat_bar));

        approx::assert_relative_eq!(
            p_dew.get::<bar>(),
            psat_bar,
            max_relative = 0.001,
        );
    }
}
