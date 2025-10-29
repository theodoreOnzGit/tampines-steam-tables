use tampines_steam_tables::interfaces::functional_programming::pt_flash_eqm;





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

    let steam_pressure = Pressure::new::<megapascal>(15.5);
    let steam_temperature = ThermodynamicTemperature::new::<degree_celsius>(300.0);
    let pri_circuit_volume = Volume::new::<cubic_meter>(239.0);
    let pri_containment_volume = Volume::new::<cubic_meter>(70_000_f64);

    // now, let us find the mass of water inside the primary circuit 
    // first by finding the specific volume of the water
    // steam quality is zero (assuming all is liquid)
    //
    let steam_quality: f64 = 0.0;

    let water_initial_specific_volume: SpecificVolume = 
        pt_flash_eqm::v_tp_eqm_two_phase(
            steam_temperature, 
            steam_pressure, 
            steam_quality);
    


}
