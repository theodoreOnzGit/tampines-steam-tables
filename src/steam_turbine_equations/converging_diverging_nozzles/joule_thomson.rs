use uom::si::f64::*;
use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;

/// we want to obtain the outlet thermodynamic state of a flow 
/// going through a joule thomson model
/// stagnation enthalpy is constant, but enthalpy will differ
///
/// In this case however, the kinetic energy is NOT negligible
pub fn get_outlet_state_joule_thomson(
    p1 : Pressure,
    h1: AvailableEnergy,
    p2: Pressure,
    mass_flowrate: MassRate,
    a1: Area,
) -> TampinesSteamTableCV {

    let ref_vol = TampinesSteamTableCV::get_ref_vol();
    let state_1 = TampinesSteamTableCV::new_from_ph(
        p1, h1, ref_vol
    );


    // so we are here going from p1 to p2 given constant mass flowrate 
    // and stagnation enthalpy

    // first, we assert that we are depressurising
    assert!(p2 < p1);
    
    let rho1 = state_1.get_rho();

    let v1: Velocity = mass_flowrate/rho1/a1;

    // stagnation enthalpy
    let h0: AvailableEnergy = h1 + 0.5 * v1 * v1;

    // Helper: Calculate mass flowrate given outlet enthalpy (p,h) flash
    fn calculate_mass_rate_at_outlet_ph(
        h0: AvailableEnergy,
        p2: Pressure,
        h2: AvailableEnergy,
        a2: Area,
    ) -> MassRate {
        // Energy equation: v₂ = √(2(h₀ - h₂))
        let v2: Velocity = (2.0 * (h0 - h2)).sqrt();
        
        // Get density from (p,h) flash
        let ref_vol = TampinesSteamTableCV::get_ref_vol();
        let state_2 = TampinesSteamTableCV::new_from_ph(p2, h2, ref_vol);
        let rho2 = state_2.get_rho();
        
        // Mass flux: G = ρv
        // Mass rate: G*a2
        let mass_rate = rho2 * v2 * a2;
        
        mass_rate
    }

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

    


    todo!()

}
/// let superheated steam flowing at 700 m/s (quite obviously supersonic)
/// 1 bar, be throttled down to 0.5 bar
///
/// temperature of steam will be 300C 
/// we shall find outlet temperature of the steam using this algorithm 
#[cfg(test)]
mod joule_thomson_test {
    use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
use crate::steam_turbine_equations::joule_thomson::get_outlet_state_joule_thomson;
    use uom::si::area::square_meter;
    use uom::si::f64::*;
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

        let state_2 = get_outlet_state_joule_thomson(
            p1, h1, p2, mass_flowrate, area
        );


    }
}

