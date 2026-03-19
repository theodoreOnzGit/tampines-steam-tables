use uom::ConstZero;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::{TampinesSteamTableCV};
use crate::prelude::functional_programming::ph_flash_eqm::v_ph_eqm;

/// basically the same code, but returns only p2 and h2
pub fn get_isentropic_nozzles_outlet_ph_point(
    p1: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    user_set_tolerance: Option<f64>,
) -> (Pressure, AvailableEnergy) {
    let (p2,h2,_rho2) = 
        get_isentropic_nozzles_outlet_ph_rho_point(
            p1, h1, mass_flowrate, a1, a2, user_set_tolerance
        );

    return (p2,h2);
}

/// given inlet and outlet areas, 
/// a1 and a2
/// the inlet conditions, p1, h1, and v1
/// mass flowrate
/// one should be able to find iteratively, p2 and h2
#[inline]
pub fn get_isentropic_nozzles_outlet_ph_rho_point(
    p1: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    user_set_tolerance: Option<f64>,
) -> (Pressure, AvailableEnergy, MassDensity) {

    let tolerance: f64;

    match user_set_tolerance {
        Some(mut user_tol) => {
            if user_tol >= 1.0 {
                user_tol = 1e-2;
            }
            tolerance = user_tol;
        },
        None => {
            tolerance = 1e-2;
        },
    }

    // now let's get the thermodynamic state 

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

    let rho1: MassDensity = state_1.get_specific_volume().recip();
    let v1: Velocity = mass_flowrate/rho1/a1;

    let mut rho2_guess: MassDensity = rho1;
    let mut p2_guess: Pressure = Pressure::ZERO;
    let mut residual = 1.0;
    let mut h2_guess = AvailableEnergy::ZERO;
    let ratio_one = Ratio::new::<ratio>(1.0);

    // now we are ready to loop 

    while residual > tolerance {

        let rho2 = rho2_guess;
        // let's define a few terms 
        //
        let rho1_a1_by_rho2_a2: Ratio = rho1 * a1 / (rho2 * a2);

        // in this part 
        // h1 - h2 = 0.5 * (rho1_a1_by_rho2_a2^2 - 1) * v1sq

        let h1_minus_h2: AvailableEnergy = 
            0.5 * (
                rho1_a1_by_rho2_a2 * rho1_a1_by_rho2_a2 
                - ratio_one
            ) * v1 * v1;

        h2_guess = h1 - h1_minus_h2;

        // then 
        // p1 a1 - p2 a2 = mass_flowrate * v1 (rho1_a1_by_rho2_a2 - 1)

        let p1a1_minus_p2a2: Force 
            = mass_flowrate * v1 * (rho1_a1_by_rho2_a2 - ratio_one);

        let p1a1: Force = p1*a1;

        let p2a2 = p1a1 - p1a1_minus_p2a2;

        p2_guess = p2a2/a2;

        // now we have new thermodynamic state

        rho2_guess = v_ph_eqm(p2_guess, h2_guess).recip();

        // then let's get residual

        residual = ((rho2 - rho2_guess)/rho2_guess).into();
        residual = residual.abs();

        let debug = false;
        if debug {
            dbg!(&(rho2_guess));
            dbg!(&(p2_guess,h2_guess));
            dbg!(&((residual,tolerance)));
        }


    };
    

    return (p2_guess, h2_guess,rho2_guess);

}

#[cfg(test)]
mod nozzles_test {
    use uom::si::f64::*;
    use uom::si::area::square_centimeter;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::bar;
    use uom::si::thermodynamic_temperature::degree_celsius;

    use crate::prelude::functional_programming::pt_flash_eqm::h_tp_eqm_single_phase;
    use crate::steam_turbine_equations::get_isentropic_nozzles_outlet_ph_rho_point;


    #[test]
    pub fn nozzles_test(){

        let a1: Area = Area::new::<square_centimeter>(5.0);
        let a2: Area = Area::new::<square_centimeter>(2.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
        let tolerance: Option<f64> = Option::None;
        // let this be steam at 10 bar 330 C
        // this is superheated steam
        let p1 = Pressure::new::<bar>(10.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(330.0);
        let h1 = h_tp_eqm_single_phase(t1, p1);
        

        let (p2,h2,rho2) = get_isentropic_nozzles_outlet_ph_rho_point(
            p1, h1, mass_flowrate, a1, a2, tolerance);

        dbg!(&(p2,h2,rho2));

        let v2: Velocity = mass_flowrate/a2/rho2;

        dbg!(&v2);

        todo!();

    }
}
