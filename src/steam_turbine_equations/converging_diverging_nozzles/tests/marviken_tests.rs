// for the Marviken tests 
//
// there was a vessel at 425 m3  (fig 3.1 on page 9)
// https://www.nrc.gov/docs/ML2005/ML20052H367.pdf
// This vessel was 

// In your test module

use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::ratio::ratio;
use uom::si::pressure::{atmosphere, kilopascal, megapascal};
use uom::si::length::millimeter;
use uom::si::velocity::meter_per_second;

use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
use crate::steam_turbine_equations::calculate_velocity_mass_flowrate_and_state_in_cd_nozzle;

/// looks like fig 8:24 on page 113 seems to be the best 
/// as most nozzles have pipe which may cause extra pressure loss and 
/// reduce flowrate.
///
/// These correspond to test 23 and 24. 
///
/// But Fig 8:24 has L/D at 0.3 
/// which makes it as pure a nozzle as it gets
///
/// In the pdf, this is on page 100 (113 of the pdf)
/// The conditions of test 23 and 24 are in table 4:2 on 
/// page 21 (page 34 of the pdf)
///
/// https://www.nrc.gov/docs/ML2005/ML20052H367.pdf
/// NUREG/CR-2671
/// MXC-301
///
/// For test 23, these are the datapoints from graphreader
/// inlet pressure (kPa), nozzle mass flux (kg/s /m2)
///
/// 3724.711,19501.04
/// 3778.902,19501.04
/// 3829.48,19501.04
/// 3887.283,19209.979
/// 3959.538,19501.04
/// 4024.566,19646.57
/// 4075.145,19792.1
/// 4125.723,19792.1
/// 4190.751,19792.1
/// 4248.555,20228.69
/// 4313.584,20374.22
/// 4385.838,20956.341
/// 4443.642,21101.871
/// 4494.22,21247.401
/// 4537.572,21392.931
/// 4580.925,22120.582
/// 4631.503,22848.233
/// 4667.63,23284.823
/// 4696.532,23721.414
/// 4747.11,24158.004
/// 4797.688,28669.439
/// 4812.139,24594.595
/// 4841.04,25322.245
/// 4877.168,31725.572
/// 4891.618,29397.089
/// 4898.844,27214.137
/// 4913.295,25467.775
/// 4942.197,32744.283
/// 4974.711,33035.343

///
/// For test 24, these are the datapoints from graphreader
/// inlet pressure (kPa), nozzle mass flux (kg/s /m2)
/// 2828.757,16735.967
/// 2868.497,17318.087
/// 2904.624,16881.497
/// 2947.977,16735.967
/// 2984.104,17172.557
/// 3027.457,22266.112
/// 3049.133,20083.16
/// 3063.584,21247.401
/// 3085.26,20519.751
/// 3121.387,30124.74
/// 3150.289,26049.896
/// 3164.74,27214.137
/// 3193.642,34636.175
/// 3215.318,33180.873
/// 3273.121,35509.356
/// 3287.572,37546.778
/// 3316.474,39147.609
/// 3352.601,40020.79
/// 3417.63,41767.152
/// 3453.757,40020.79
/// 3576.59,43513.514
/// 3612.717,44386.694
/// 3634.393,43513.514
/// 3706.647,44386.694
/// 3778.902,45405.405
/// 3822.254,44823.285
/// 3901.734,45405.405
/// 3916.185,46424.116
/// 3959.538,46133.056
/// 4060.694,47442.827
/// 4075.145,48898.129
/// 4161.85,51808.732
/// 4255.78,50935.551
/// 4306.358,51808.732
/// 4356.936,52827.443
/// 4421.965,53700.624
/// 4515.896,53409.563
/// 4580.925,54573.805
/// 4703.757,54137.214
/// 4772.399,56611.227
///
///
///
#[test]
#[ignore="skip first, Marviken is more complex"]
fn validate_against_marviken_test_24() {


    // let's first have our datapoints 

    let pressure_and_mass_flux_vec_test_24 = vec![
        (2828.757,16735.967),
        (2868.497,17318.087),
        (2904.624,16881.497),
        (2947.977,16735.967),
        (2984.104,17172.557),
        (3027.457,22266.112),
        (3049.133,20083.16),
        (3063.584,21247.401),
        (3085.26,20519.751),
        (3121.387,30124.74),
        (3150.289,26049.896),
        (3164.74,27214.137),
        (3193.642,34636.175),
        (3215.318,33180.873),
        (3273.121,35509.356),
        (3287.572,37546.778),
        (3316.474,39147.609),
        (3352.601,40020.79),
        (3417.63,41767.152),
        (3453.757,40020.79),
        (3576.59,43513.514),
        (3612.717,44386.694),
        (3634.393,43513.514),
        (3706.647,44386.694),
        (3778.902,45405.405),
        (3822.254,44823.285),
        (3901.734,45405.405),
        (3916.185,46424.116),
        (3959.538,46133.056),
        (4060.694,47442.827),
        (4075.145,48898.129),
        (4161.85,51808.732),
        (4255.78,50935.551),
        (4306.358,51808.732),
        (4356.936,52827.443),
        (4421.965,53700.624),
        (4515.896,53409.563),
        (4580.925,54573.805),
        (4703.757,54137.214),
        (4772.399,56611.227), 
        ];

    // --- Step 2: Define Geometry from Figure 3.2.2 ---
    let nozzle_diameter = Length::new::<millimeter>(500.0);
    let nozzle_area = std::f64::consts::PI * (nozzle_diameter * nozzle_diameter / 4.0);
    let a_throat = nozzle_area;
    let a_exit =  nozzle_area; // It's a converging nozzle, so throat area = exit area
    // The back pressure is atmospheric, as they are venting to a large containment vessel.
    // The back pressure is atmospheric, as they are venting to a large containment vessel.

    for (pressure_kpa_ptr, mass_flux_kg_per_s_m2) 
        in pressure_and_mass_flux_vec_test_24.iter() {

        let initial_vessel_pressure = Pressure::new::<megapascal>(4.95);
        // --- Step 3: Get Initial State from your Steam Tables ---
        let ref_vol = TampinesSteamTableCV::get_ref_vol();
        let p2 = Pressure::new::<atmosphere>(1.0);
        let p1 = Pressure::new::<kilopascal>(*pressure_kpa_ptr);
        let t1 = TampinesSteamTableCV::try_get_tsat(p1).unwrap() - 
            TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(33.0);
        let state_1 = TampinesSteamTableCV::new_from_tp_quality_1(t1, p1, ref_vol);
        let h1 = state_1.get_specific_enthalpy();

        // The velocity inside the huge Marviken pressure vessel is effectively zero.
        let v1 = Velocity::new::<meter_per_second>(0.0);

        // --- Step 4: Call Your Master Function ---
        let (v_out, m_dot_out, state_out) = 
            calculate_velocity_mass_flowrate_and_state_in_cd_nozzle(
                p1, h1, v1, a_throat, a_exit, p2
            );

        
        // --- Step 5: Compare to the Experimental Result from Table 5.1.1 ---
        let experimental_mass_flux = 
            MassFlux::new::<kilogram_per_square_meter_second>(
                *mass_flux_kg_per_s_m2
            );

        let calculated_mass_flux = 
            m_dot_out/a_throat;

        //dbg!(&(initial_vessel_pressure,p1,p2));
            

        println!("Calculated Mass Flux: {:.2} kg/(m2 s)", calculated_mass_flux.get::<kilogram_per_square_meter_second>());
        println!("Experimental Mass Flux: {:.2} kg/(m2 s)", experimental_mass_flux.get::<kilogram_per_square_meter_second>());

        // Use a relative difference assertion. Don't expect a perfect match!
        // Getting within 10-15% would be a fantastic result for a 1D model.
        let relative_difference = ((calculated_mass_flux - experimental_mass_flux) / experimental_mass_flux).abs();

        println!("Relative Difference: {:.2}%", relative_difference.get::<ratio>() * 100.0);

        // Assert that your result is within a reasonable tolerance, e.g., 20%
        //assert!(relative_difference.get::<ratio>() < 0.20);

    }

    todo!()

    



}
