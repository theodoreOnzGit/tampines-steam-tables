








// this example is meant to show some capabilities of tampines-steam-tables
//
// suppose we have a pressurised water reactor, 15.5 MPa, 300 degrees C 
//
// sudden depressurisation occurs into a containment vessel
//
// Malet, J., Parduba, Z., Mimouni, S., & Travis, J. (2014). Achievements 
// of spray activities in nuclear reactor containments during the 
// last decade. Annals of Nuclear Energy, 74, 134-142.
//
// These containments have on the order of 70,000 m3 of volume
//
// For the reactor pressure vessel
//
// https://www-pub.iaea.org/MTCD/Publications/PDF/TCS-22_web.pdf
// No, T. C. S. (2003). Pressurized water reactor simulator.
//
// typical primary circuit volume is about 239 m3,
fn main(){


    println!("starting depressurisation test");

    // let's start with the state of steam first
    // note that all quantities require the use of units of measure (uom)
    // eg. temperature and pressure, 
    // these must be dimensioned

    use uom::si::pressure::megapascal;
    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::volume::cubic_meter;
    use uom::si::specific_volume::cubic_meter_per_kilogram;
    use uom::si::pressure::bar;
    use uom::si::ratio::ratio;

    let steam_pressure = Pressure::new::<megapascal>(15.5);
    let steam_temperature = ThermodynamicTemperature::new::<degree_celsius>(300.0);
    let pri_circuit_volume = Volume::new::<cubic_meter>(239.0);
    let pri_containment_volume = Volume::new::<cubic_meter>(70_000_f64);

    // now, let us find the mass of water inside the primary circuit 
    // first by finding the specific volume of the water
    // steam quality is zero (assuming all is liquid)
    //
    let steam_quality: f64 = 0.0;

    // this is in m3/kg
    // to obtain the specific volume, we use the steam temperature,
    // pressure and quality to obtain it, 
    // using a tp(temperature and pressure) equilibrium calculation
    // ie, waiting for the steam to reach same temperature and pressure 
    // throughout, ie. thermodynamic equilibrium
    use tampines_steam_tables::interfaces::functional_programming::pt_flash_eqm;
    let water_initial_specific_volume: SpecificVolume = 
        pt_flash_eqm::v_tp_eqm_two_phase(
            steam_temperature, 
            steam_pressure, 
            steam_quality);
    // we can find the mass of water inside the RPV 
    let water_mass: Mass = 
        pri_circuit_volume/water_initial_specific_volume;

    // given the mass of water, we can now
    // find the specific volume at the end

    let water_final_specific_volume_estimate: SpecificVolume = 
        pri_containment_volume/water_mass;

    // the specific volume is around 0.403 m3/kg
    assert_eq!(water_final_specific_volume_estimate,
        SpecificVolume::new::<cubic_meter_per_kilogram>(0.40314065396255444)
    );

    // now what is the final pressure of this system?
    // we don't have much information except for specific volume,
    // 
    //
    // Also, we know the enthalpy of the water at the beginning should 
    // equal the enthalpy of the water at the end (if we don't factor 
    // in decay heat or reactor heat)
    //
    // note AvailableEnergy is in Energy/mass

    use tampines_steam_tables::interfaces::functional_programming::pt_flash_eqm::h_tp_eqm_two_phase;
    let water_enthalpy: AvailableEnergy
        = h_tp_eqm_two_phase(
            steam_temperature, 
            steam_pressure, 
            steam_quality
        );

    // the enthalpy at the end is the same,
    // however, there is depressurisation
    // to what degree, we don't know 
    // however, the more depressurisation there is, the less the 
    // specific volume 
    //
    // there is no direct method to calculate this, however, we can  
    // do so iteratively.
    //
    // this means we guess and check the pressure iteratively until 
    // we obtain the correct pressure
    //
    // Now there are pre-coded Rust methods for this, 
    // for other crates
    // eg. bisection
    //
    // but I'm going to demonstrate it here

    
    // suppose we have a guessed pressure of 2 bar
    let two_bar_guessed_pressure = Pressure::new::<bar>(2.0);

    // let's calculate a specific volume value:
    // we shall use the ph, flash algorithm, as in pressure and enthalpy 
    // equilibrium calculation
    use tampines_steam_tables::interfaces::functional_programming::ph_flash_eqm::v_ph_eqm;
    let two_bar_specific_volume: SpecificVolume 
        = v_ph_eqm(two_bar_guessed_pressure, water_enthalpy);

    println!("the specific volume at 2 bar is:");
    println!("{:?}",two_bar_specific_volume);
    assert_eq!(two_bar_specific_volume,
        SpecificVolume::new::<cubic_meter_per_kilogram>(0.33577033964413827)
    );

    // note that this specific volume is lower than the final specific volume (assuming the steam fills the whole containment)

    let upper_bound_pressure = two_bar_guessed_pressure;

    // now for larger specific volume, we lower the pressure 
    // say 0.5 bar
    //
    let half_bar_guessed_pressure = Pressure::new::<bar>(0.5);

    let half_bar_specific_volume: SpecificVolume 
        = v_ph_eqm(half_bar_guessed_pressure, water_enthalpy);

    println!("the specific volume at 0.5 bar is:");
    println!("{:?}",half_bar_specific_volume);
    assert_eq!(half_bar_specific_volume,
        SpecificVolume::new::<cubic_meter_per_kilogram>(1.4024427855694763)
    );

    // so this specific volume is bigger than the 0.403 m3/kg roughly 
    // so the actual pressure is somewhere between 0.5 to 2 bar
    let lower_bound_pressure = half_bar_guessed_pressure;

    // to do the solution, we can do a bisection iteration algorithm 
    // like so:

    fn bisection_algo(lower_bound_pressure: Pressure,
        upper_bound_pressure: Pressure,
        water_enthalpy: AvailableEnergy,
        target_specific_volume: SpecificVolume) -> (Pressure, SpecificVolume) {

        // first we specify a tolerance 

        let tolerance = Ratio::new::<ratio>(1e-5);

        // let's have some initial guesses

        let mut guessed_pressure = 0.5 * (upper_bound_pressure +
            lower_bound_pressure);
        let mut upper_pressure_bracket: Pressure = upper_bound_pressure;
        let mut lower_pressure_bracket: Pressure = lower_bound_pressure;

        let mut guess_specific_volume 
            = v_ph_eqm(guessed_pressure, water_enthalpy);

        let mut error: Ratio = (guess_specific_volume - 
            target_specific_volume)/target_specific_volume;

        while error.abs() > tolerance {

            // first check if the specific volume is more or less than 
            // the target specific volume 

            if guess_specific_volume > target_specific_volume {
                // in the case it is greater than the target specific 
                // volume 
                // that means pressure is too low 
                lower_pressure_bracket = guessed_pressure;

            } else {
                // in the case it is lower than the guessed specific 
                // volume, that means pressure is to high 
                upper_pressure_bracket = guessed_pressure;

            }

            // we'll try to guess the pressure again 
            guessed_pressure = 0.5 * (upper_pressure_bracket 
                + lower_pressure_bracket);
            // and then the specific volume
            guess_specific_volume 
                = v_ph_eqm(guessed_pressure, water_enthalpy);
            // recalculate the error
            error = (guess_specific_volume - 
                target_specific_volume)/target_specific_volume;
            // so if error is sufficiently low, the guessed pressure 
            // will be the correct pressure
            // we exit the loop

        };

        return (guessed_pressure, guess_specific_volume);

    }


    let target_specific_volume = water_final_specific_volume_estimate;
    let (final_pressure, final_specific_volume) = bisection_algo(
        lower_bound_pressure, 
        upper_bound_pressure, 
        water_enthalpy, 
        target_specific_volume);


    println!("the final pressure is:");
    println!("{:?}",final_pressure);
    

    println!("the final specific volume is:");
    println!("{:?}",final_specific_volume);



    // hence, after depressurisation, we arrive at a final pressure of 
    // about 1.68 bar
    // what is the water quality?
    // as in how much liquid water remains?
    // we can use the quality calculation

    use tampines_steam_tables::interfaces::functional_programming::ph_flash_eqm::x_ph_flash;
    let final_quality = x_ph_flash(final_pressure, water_enthalpy);

    println!("the final quality is:");
    println!("{:?}",final_quality);
    assert_eq!(final_quality,
        0.38614681771539744
    );
    // which means we have a saturated vapour/liquid equilibrium
    // quality is (mass of vap)/(mass of vap + mass of liq)
    // note that mass of vapour and liquid is the initial mass of 
    // pressurised water
    //
    let liquid_mass: Mass = 
        final_quality * water_mass;

    println!("the final mass of liquid water is:");
    println!("{:?}",liquid_mass);
    println!("out of an initial water coolant inventory of");
    println!("{:?}",water_mass);
    // we have about 67,000 kg of water remaining, after the depressurisation

}
