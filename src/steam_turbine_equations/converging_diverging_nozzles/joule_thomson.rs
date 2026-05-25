use uom::si::{f64::*, ratio::ratio};
use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;

/// we want to obtain the outlet thermodynamic state of a flow 
/// going through a joule thomson model
/// stagnation enthalpy is constant, but enthalpy will differ
///
/// In this case however, the kinetic energy is NOT negligible
pub fn get_outlet_velocity_and_state_joule_thomson(
    p1 : Pressure,
    h1: AvailableEnergy,
    p2: Pressure,
    mass_flowrate_ref: MassRate,
    a1: Area,
) -> (Velocity, TampinesSteamTableCV) {

    let ref_vol = TampinesSteamTableCV::get_ref_vol();
    let state_1 = TampinesSteamTableCV::new_from_ph(
        p1, h1, ref_vol
    );


    // so we are here going from p1 to p2 given constant mass flowrate 
    // and stagnation enthalpy

    // first, we assert that we are depressurising
    assert!(p2 < p1);
    
    let rho1 = state_1.get_rho();

    let v1: Velocity = mass_flowrate_ref/rho1/a1;

    // stagnation enthalpy
    let h0: AvailableEnergy = h1 + 0.5 * v1 * v1;



    // Helper: calculate mass flowrate using outlet enthalpy (p,h) flash 
    // using velocit as input
    fn calculate_mass_rate_and_state_at_outlet_ph_velocity(
        h0: AvailableEnergy,
        p2: Pressure,
        v2: Velocity,
        a2: Area,
    ) -> (MassRate, TampinesSteamTableCV) {
        // Energy equation: v₂ = √(2(h₀ - h₂))
        //
        // we use: h2 = h0 - 0.5 * v2^2
        let h2: AvailableEnergy = h0 - 0.5 * v2 * v2;
        
        // Get density from (p,h) flash
        let ref_vol = TampinesSteamTableCV::get_ref_vol();
        let state_2 = TampinesSteamTableCV::new_from_ph(p2, h2, ref_vol);
        let rho2 = state_2.get_rho();
        
        // Mass flux: G = ρv
        // Mass rate: G*a2
        let mass_rate = rho2 * v2 * a2;
        
        (mass_rate, state_2)
    }

    // For such a case, in expansion throttling,
    // we expect that density decreases (specific volume increases)
    // so that 
    // velocity increases under constant area
    // from this, we know that there is an upper limit to what densities 
    // are available
    //
    // we will start at 5% density reduction at each step, until the 
    // sign changes

    let mut test_rho = rho1;
    let mut mass_rate_error_initial = 0.0;
    let mut v_upper_limit: Velocity = v1;
    let mut v_lower_limit: Velocity = v1;

    let max_iterations = 50;
    let debug = false;
    let mut relative_error: f64 = 1.0;
    const TOLERANCE: f64 = 0.0001;  // 0.01% tolerance


    for i in 0..max_iterations {

        // first let's test the mass flowrate and state 
        // by obtaining the velocity

        let v2: Velocity = mass_flowrate_ref/test_rho/a1;
        let a2 = a1;

        let (mass_rate_test, state_2_test) = 
            calculate_mass_rate_and_state_at_outlet_ph_velocity(
                h0, p2, v2, a2
            );
        
        // check residual 
        //
        relative_error = 
            ((mass_rate_test - mass_flowrate_ref)/mass_flowrate_ref)
            .get::<ratio>();

        // for ONLY the first iteration, we check the initial mass rate 
        // error 

        if i == 0 {
            mass_rate_error_initial = relative_error;
        }

        

        // mutate the test rho, to reduce it

        if debug {
            dbg!(&(test_rho,relative_error,state_2_test));
        }

        // if there is sign change, break out 
        // we capture the upper bound velocity

        if relative_error * mass_rate_error_initial <= 0.0 {
            v_upper_limit = v2;
            break;
        }
        // if there is no sign change, update the lower bound velocity 

        v_lower_limit = v2;
        test_rho *= 0.90;

    }
    // now for sure, the velocity lies between the upper and lower bound
    // we can use regula falsi to finish this job

    if debug {
        dbg!(&(v_lower_limit, v_upper_limit));
    }

    let root_finder = |v_test: Velocity| -> (MassRate, TampinesSteamTableCV) {

        let (mass_flowrate_calc, outlet_state_with_shocks) = 
            calculate_mass_rate_and_state_at_outlet_ph_velocity(
                h0, p2, v_test, a1);

        let error = mass_flowrate_calc - mass_flowrate_ref;


        (error, outlet_state_with_shocks)
    };
    

    let (mut error_lower_limit, _state_lower_limit) = 
        root_finder(v_lower_limit);
    let (mut error_upper_limit, _state_upper_limit) = 
        root_finder(v_upper_limit);

    if error_lower_limit.value * error_upper_limit.value >= 0.0 {
        panic!("bounds are same sign!");
    }


    let mut iteration = 0;
    let mut v_test;

    // this is regula falsi
    while relative_error.abs() > TOLERANCE && iteration < max_iterations {

        // using secant formula
        v_test = 
            v_upper_limit - (v_upper_limit - v_lower_limit) * 
            error_upper_limit/(error_upper_limit - error_lower_limit);

        // check mass flowrate using velocity
        let (test_error, test_outlet_state) = root_finder(v_test);
        // update relative error
        relative_error = (test_error/mass_flowrate_ref).get::<ratio>();

        if relative_error.abs() < TOLERANCE {
            return (v_test, test_outlet_state);
        }

        // update bounds, keep root bracketed 

        if error_lower_limit.value * test_error.value < 0.0 {

            v_upper_limit = v_test;
            error_upper_limit = test_error;
        } else {

            v_lower_limit = v_test;
            error_lower_limit = test_error;
        }


        if debug {
            dbg!(&(v_lower_limit,v_upper_limit));
        }
        iteration += 1;

    }

    panic!("Joule-Thomson algorithm: failed to converge");

}
/// let superheated steam flowing at 700 m/s (quite obviously supersonic)
/// 1 bar, be throttled down to 0.5 bar
///
/// temperature of steam will be 300C 
/// we shall find outlet temperature of the steam using this algorithm 
#[cfg(test)]
mod joule_thomson_test {
    use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
    use crate::steam_turbine_equations::joule_thomson::get_outlet_velocity_and_state_joule_thomson;
    use uom::si::area::square_meter;
    use uom::si::f64::*;
    use uom::si::ratio::ratio;
    use uom::si::velocity::meter_per_second;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::pressure::atmosphere;


    #[test]
    fn superheated_steam(){
        let ref_vol = TampinesSteamTableCV::get_ref_vol();
        let v1 = Velocity::new::<meter_per_second>(700.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(300.0);
        let p1 = Pressure::new::<atmosphere>(1.0);
        let p2 = Pressure::new::<atmosphere>(0.5);
        let area = Area::new::<square_meter>(0.5);

        let state_1 = TampinesSteamTableCV::new_from_tp_quality_1(
            t1, p1, ref_vol
        );

        let mass_flowrate: MassRate = 
            state_1.get_rho() * v1 * area;
        let h1 = state_1.get_specific_enthalpy();

        let (v2, state_2) = get_outlet_velocity_and_state_joule_thomson(
            p1, h1, p2, mass_flowrate, area
        );

        assert!(v2 > v1);

        dbg!(&(v1,state_1));
        dbg!(&(v2,state_2));
        let h2 = state_2.get_specific_enthalpy();

        let h0_ref = 0.5 * v1 * v1 + h1;
        let h0_test = 0.5 * v2 * v2 + h2;

        let enthalpy_error = ((h0_test - h0_ref)/h0_ref).get::<ratio>();

        assert!(enthalpy_error < 1e-5);

        // we also just check v2 
        // it is higher than 700 m/s, this just gives an idea of the 
        // acceleration involved
        //
        // and is a regression test

        approx::assert_relative_eq!(
            v2.get::<meter_per_second>(),
            1035.67,
            max_relative=1e-4
        );

    }
}

