use uom::si::{f64::*, ratio::ratio, volume::cubic_meter};

use crate::prelude::{TampinesSteamTableCV, functional_programming::ph_flash_eqm::v_ph_eqm};

/// given inlet and outlet areas, 
/// a1 and a2
/// the inlet conditions, p1, h1, and v1
/// mass flowrate
/// one should be able to find iteratively, p2 and h2
pub fn get_isentropic_nozzles_outlet_ph_point(
    p1: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    v1: Velocity,
    tolerance: Option<f64>,
) -> (Pressure, AvailableEnergy) {

    let tol: f64;

    match tolerance {
        Some(mut user_tol) => {
            if user_tol >= 1.0 {
                user_tol = 1e-2;
            }
            tol = user_tol;
        },
        None => {
            tol = 1e-2;
        },
    }

    // now let's get the thermodynamic state 

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

    let rho1: MassDensity = state_1.get_specific_volume().recip();

    let mut rho2_guess: MassDensity = rho1;
    let mut p2_guess: Pressure;
    let mut residual = 1.0;
    let mut h2_guess;
    let ratio_one = Ratio::new::<ratio>(1.0);

    // now we are ready to loop 

    while residual > tol {

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

        let debug = true;
        if debug {
            dbg!(&(rho2_guess));
            dbg!(&(p2_guess,h2_guess));
        }


    };
    

    return (p2_guess, h2_guess);

}
