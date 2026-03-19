use uom::ConstZero;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::hs_flash_eqm::v_hs_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::v_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::prelude::{TampinesSteamTableCV};

/// basically the same code, but returns only p2 and h2
///
/// Note: this uses (h,s) equations, so it can only cover part of the steam 
/// table
pub fn get_isentropic_nozzles_outlet_ph_point(
    p1: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    user_set_tolerance: Option<f64>,
) -> (Pressure, AvailableEnergy) {
    let (p2,h2,_rho2) = 
        get_isentropic_nozzles_outlet_ph_rho_point_hs_algo(
            p1, h1, mass_flowrate, a1, a2, user_set_tolerance
        );

    return (p2,h2);
}

/// given inlet and outlet areas, 
/// a1 and a2
/// the inlet conditions, p1, h1, and v1
/// mass flowrate
/// one should be able to find iteratively, p2 and h2
///
/// these are very simple iterative equations that back-substitute the 
/// density 
///
/// Note: this uses (h,s) equations, so it can only cover part of the steam 
/// table
#[inline]
pub fn get_isentropic_nozzles_outlet_ph_rho_point_hs_algo(
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
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();
    // note: this is isentropic, hence s2 = s1;
    let s2 = s1;
    let v1: Velocity = mass_flowrate/rho1/a1;

    let mut rho2_guess: MassDensity = rho1;
    let mut p2_guess: Pressure = Pressure::ZERO;
    let mut rho_residual = 1.0;
    let mut h2_guess = AvailableEnergy::ZERO;
    let ratio_one = Ratio::new::<ratio>(1.0);
    let p1a1: Force = p1*a1;

    // now we are ready to loop 

    while rho_residual > tolerance {

        let rho2 = rho2_guess;
        // let's define a few terms 
        //
        let mut rho1_a1_by_rho2_a2: Ratio = rho1 * a1 / (rho2 * a2);

        // in this part 
        // h1 - h2 = 0.5 * (rho1_a1_by_rho2_a2^2 - 1) * v1sq

        let h1_minus_h2: AvailableEnergy = 
            0.5 * (
                rho1_a1_by_rho2_a2 * rho1_a1_by_rho2_a2 
                - ratio_one
            ) * v1 * v1;

        h2_guess = h1 - h1_minus_h2;

        rho2_guess = v_hs_eqm(h2_guess, s2).recip();

        rho1_a1_by_rho2_a2 = rho1 * a1 / (rho2_guess * a2);


        // then 
        // p1 a1 - p2 a2 = mass_flowrate * v1 (rho1_a1_by_rho2_a2 - 1)

        let p1a1_minus_p2a2: Force 
            = mass_flowrate * v1 * (rho1_a1_by_rho2_a2 - ratio_one);


        let p2a2 = p1a1 - p1a1_minus_p2a2;

        p2_guess = p2a2/a2;

        // now we have new thermodynamic state

        //rho2_guess = v_ph_eqm(p2_guess, h2_guess).recip();

        // then let's get residual

        rho_residual = ((rho2 - rho2_guess)/rho2_guess).into();
        rho_residual = rho_residual.abs();

        let debug = false;
        if debug {
            dbg!(&(rho2_guess));
            dbg!(&(p2_guess,h2_guess));
            dbg!(&((rho_residual,tolerance)));
        }


    };
    

    return (p2_guess, h2_guess,rho2_guess);

}


/// given inlet and outlet areas, 
/// a1 and a2
/// the inlet conditions, p1, h1, and v1
/// mass flowrate
/// one should be able to find iteratively, p2 and h2
///
/// these are very simple iterative equations that back-substitute the 
/// density 
///
/// Note: this uses (p,s) equations, so it 
/// can cover all of the steam table
/// table
#[inline]
pub fn get_isentropic_nozzles_outlet_ph_rho_point_ps_algo(
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
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();
    // note: this is isentropic, hence s2 = s1;
    let s2 = s1;
    let v1: Velocity = mass_flowrate/rho1/a1;

    let p1a1: Force = p1*a1;
    let mut rho2_guess: MassDensity = rho1;
    let mut p2_guess: Pressure = Pressure::ZERO;
    let mut rho_residual = 1.0;
    let mut h2_guess = AvailableEnergy::ZERO;
    let ratio_one = Ratio::new::<ratio>(1.0);

    // now we are ready to loop 

    while rho_residual > tolerance {

        let rho2 = rho2_guess;
        // let's define a few terms 
        //
        let mut rho1_a1_by_rho2_a2: Ratio = rho1 * a1 / (rho2 * a2);
        // in the (p,s) algorithm, we first guess p2
        // p1 a1 - p2 a2 = mass_flowrate * v1 (rho1_a1_by_rho2_a2 - 1)

        let mut p1a1_minus_p2a2: Force 
            = mass_flowrate * v1 * (rho1_a1_by_rho2_a2 - ratio_one);


        let mut p2a2 = p1a1 - p1a1_minus_p2a2;

        p2_guess = p2a2/a2;

        rho2_guess = v_ps_eqm(p2_guess, s2).recip();
        rho1_a1_by_rho2_a2 = rho1 * a1 / (rho2_guess * a2);

        // with this p2_guess, we can guess the density at the new state
        // in this part 
        // h1 - h2 = 0.5 * (rho1_a1_by_rho2_a2^2 - 1) * v1sq

        let h1_minus_h2: AvailableEnergy = 
            0.5 * (
                rho1_a1_by_rho2_a2 
                + ratio_one
            ) * v1 * p1a1_minus_p2a2/mass_flowrate
            ;

        h2_guess = h1 - h1_minus_h2;


        // need to update p2_guess after

        p1a1_minus_p2a2
            = mass_flowrate * v1 * (rho1_a1_by_rho2_a2 - ratio_one);


        p2a2 = p1a1 - p1a1_minus_p2a2;

        p2_guess = p2a2/a2;


        // now we have new thermodynamic state

        //rho2_guess = v_ph_eqm(p2_guess, h2_guess).recip();

        // then let's get residual

        rho_residual = ((rho2 - rho2_guess)/rho2_guess).into();
        rho_residual = rho_residual.abs();

        let debug = true;
        if debug {
            dbg!(&(rho2_guess));
            dbg!(&(p2_guess,h2_guess));
            dbg!(&((rho_residual,tolerance)));
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
    use uom::si::ratio::ratio;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::volume::cubic_meter;

    use crate::prelude::TampinesSteamTableCV;
    use crate::prelude::functional_programming::pt_flash_eqm::h_tp_eqm_single_phase;
    use crate::steam_turbine_equations::{get_isentropic_nozzles_outlet_ph_rho_point_hs_algo, get_isentropic_nozzles_outlet_ph_rho_point_ps_algo};


    #[test]
    pub fn nozzles_test_hs_algo(){

        let a1: Area = Area::new::<square_centimeter>(5.0);
        let a2: Area = Area::new::<square_centimeter>(4.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
        let tolerance: Option<f64> = Option::None;
        // let this be steam at 10 bar 330 C
        // this is superheated steam
        let p1 = Pressure::new::<bar>(10.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(330.0);
        let h1 = h_tp_eqm_single_phase(t1, p1);
        

        let (p2,h2,rho2) = get_isentropic_nozzles_outlet_ph_rho_point_hs_algo(
            p1, h1, mass_flowrate, a1, a2, tolerance);


        let v2: Velocity = mass_flowrate/a2/rho2;

        dbg!(&(p2,h2,rho2,v2));

        let m2: MassRate = rho2*v2*a2;

        let mass_flowrate_residual: f64 = 
            (1.0 - (m2/mass_flowrate).get::<ratio>()).abs();

        // to test for internal consistency, we do mass, energy and 
        // momentum balance before and after the nozzle

        // for conservation of mass
        // assert that the mass flowrate is within 0.1%
        assert!(mass_flowrate_residual < 1e-3);

        // for conservation of energy 
        // assert that stagnation enthalpy is within 0.1%

        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_1: TampinesSteamTableCV = 
            TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

        let rho1: MassDensity = state_1.get_specific_volume().recip();
        let v1 = mass_flowrate/rho1/a1;
        let stagnation_enthalpy_1: AvailableEnergy = 
            h1 + 0.5 * v1*v1;

        let stagnation_enthalpy_2: AvailableEnergy = 
            h2 + 0.5 *v2*v2;

        let stagnation_enthalpy_residual: f64 = 
            (1.0 - (stagnation_enthalpy_1/stagnation_enthalpy_2) 
            .get::<ratio>())
            .abs();


        assert!(stagnation_enthalpy_residual < 1e-3);

        // thirdly, we must assert the force balance for conservation 
        // of momentum
        // that is 
        // p1a1 - p2a2 = m_flowrate (v2 -v1)

        let nozzle_force: Force = p1*a1 - p2*a2;

        let momentum_change_rate: Force = mass_flowrate *(v2-v1);

        let force_bal_residual: f64 = 
            (1.0 - (nozzle_force/momentum_change_rate).get::<ratio>()).abs();

        dbg!(&(force_bal_residual,mass_flowrate_residual,stagnation_enthalpy_residual));

        assert!(force_bal_residual < 1e-3);
        todo!();

    }
    #[test]
    pub fn nozzles_test_ps_algo(){

        let a1: Area = Area::new::<square_centimeter>(5.0);
        let a2: Area = Area::new::<square_centimeter>(4.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
        let tolerance: Option<f64> = Option::None;
        // let this be steam at 10 bar 330 C
        // this is superheated steam
        let p1 = Pressure::new::<bar>(10.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(330.0);
        let h1 = h_tp_eqm_single_phase(t1, p1);
        

        let (p2,h2,rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo(
            p1, h1, mass_flowrate, a1, a2, tolerance);


        let v2: Velocity = mass_flowrate/a2/rho2;

        dbg!(&(p2,h2,rho2,v2));

        let m2: MassRate = rho2*v2*a2;

        let mass_flowrate_residual: f64 = 
            (1.0 - (m2/mass_flowrate).get::<ratio>()).abs();

        // to test for internal consistency, we do mass, energy and 
        // momentum balance before and after the nozzle

        // for conservation of mass
        // assert that the mass flowrate is within 0.1%
        assert!(mass_flowrate_residual < 1e-3);

        // for conservation of energy 
        // assert that stagnation enthalpy is within 0.1%

        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_1: TampinesSteamTableCV = 
            TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

        let rho1: MassDensity = state_1.get_specific_volume().recip();
        let v1 = mass_flowrate/rho1/a1;
        let stagnation_enthalpy_1: AvailableEnergy = 
            h1 + 0.5 * v1*v1;

        let stagnation_enthalpy_2: AvailableEnergy = 
            h2 + 0.5 *v2*v2;

        let stagnation_enthalpy_residual: f64 = 
            (1.0 - (stagnation_enthalpy_1/stagnation_enthalpy_2) 
            .get::<ratio>())
            .abs();


        assert!(stagnation_enthalpy_residual < 1e-3);

        // thirdly, we must assert the force balance for conservation 
        // of momentum
        // that is 
        // p1a1 - p2a2 = m_flowrate (v2 -v1)

        let nozzle_force: Force = p1*a1 - p2*a2;

        let momentum_change_rate: Force = mass_flowrate *(v2-v1);

        let force_bal_residual: f64 = 
            (1.0 - (nozzle_force/momentum_change_rate).get::<ratio>()).abs();

        dbg!(&(force_bal_residual,mass_flowrate_residual,stagnation_enthalpy_residual));

        assert!(force_bal_residual < 1e-3);

        todo!();

    }
}
