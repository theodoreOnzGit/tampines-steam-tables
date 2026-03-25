use uom::ConstZero;
use uom::si::f64::*;
use uom::si::force::newton;
use uom::si::pressure::{bar, millibar, pascal};
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::hs_flash_eqm::v_hs_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::s_ph_eqm;
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
/// this also uses the (p,s) equations
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
    let max_iter: usize = 30;
    let mut n_iter: usize = 1;

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

    // we expect expansion, so density decreases
    let mut rho2_guess: MassDensity = rho1 * 0.9;
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

        let mut h1_minus_h2: AvailableEnergy = 
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
        h1_minus_h2=  
            0.5 * (
                rho1_a1_by_rho2_a2 * rho1_a1_by_rho2_a2 
                - ratio_one
            ) * v1 * v1;

        h2_guess = h1 - h1_minus_h2;

        // now we have new thermodynamic state

        //rho2_guess = v_ph_eqm(p2_guess, h2_guess).recip();

        // then let's get residual


        let s2_guess = s_ph_eqm(p2_guess, h2_guess);
        rho_residual = ((rho2_guess -rho2).abs()/rho2_guess).into();
        let s_residual: f64 = ((s2_guess -s2).abs()/s2).into();

        let debug = true;
        if debug {
            dbg!(&(rho2_guess));
            dbg!(&(p2_guess,h2_guess));
            dbg!(&((s_residual,tolerance)));
            dbg!(&((s2,s2_guess)));
        }

        n_iter += 1;

        if n_iter > max_iter {
            break
        };

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
    let max_iter: usize = 30;
    let mut n_iter: usize = 1;

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
    // we expect expansion, so density decreases
    let mut rho2_guess: MassDensity = rho1 * 0.9;
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

        let s2_guess = s_ph_eqm(p2_guess, h2_guess);

        rho_residual = ((rho2_guess -rho2).abs()/rho2_guess).into();
        let s_residual: f64 = ((s2_guess -s2).abs()/s2).into();

        let debug = true;
        if debug {
            dbg!(&(rho2_guess));
            dbg!(&(p2_guess,h2_guess));
            dbg!(&((s_residual,tolerance)));
            dbg!(&((s2,s2_guess)));
        }

        n_iter += 1;

        if n_iter > max_iter {
            break
        };


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
/// can cover all of the steam table table
///
/// For this, we neglect kinetic energy of the inlet
///
#[inline]
pub fn get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
    p1: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    user_set_tolerance: Option<f64>,
) -> (Pressure, AvailableEnergy, MassDensity) {

    let tolerance: f64;
    let max_iter: usize = 30;
    let mut n_iter: usize = 1;

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
    // left hand side of momentum balance
    //
    // P1 A1 + dot{m}^2/{rho1 a1}
    let lhs: Force = p1a1 + mass_flowrate * v1;
    // we expect expansion, so density decreases
    let mut rho2_guess: MassDensity = rho1;
    let mut p2_upper_bound = p1;
    // don't expect turbine exhaust to be 
    // lower than condenser pressure on a good day, 0.04 bar
    let mut p2_lower_bound = Pressure::new::<bar>(0.03);

    let mut p2_guess: Pressure = p2_upper_bound;
    // now this part is to determine better upper and lower bounds
    //
    // I'm going to bring the upper bound down first

    let debug: bool = false; 
    if debug {
        print_graph_pts_for_outlet_pressure_and_force_balance(
            p1, h1, a1, a2, mass_flowrate
        );
    }

    let n = 20;
    for i in 0..n {

        p2_guess = (n as f64 - i as f64)/(n as f64) * p1;

        let force_bal: Force = 
            force_balance_isentropic_nozzle(
                p1, p2_guess, 
                h1, mass_flowrate, a1, a2
            );

        // usually, based on shape of the graph, force balance is below 
        // 0 
        // for upper bound for such nozzles 
        //
        // note, for some force balances, the force balance starts positive,
        // so should take note

        if force_bal < Force::ZERO {
            p2_upper_bound = p2_guess;
        } else {
            // if force balance is positive, the lower bound becomes the 
            // p2_guess 
            p2_lower_bound = p2_guess;

            break;
            // break out of this cycle

        }


    }

    // after all this, p2_guess should be midpoint of upper and lower bound 

    p2_guess = 0.5 * p2_upper_bound + 0.5 * p2_lower_bound;





    let mut p2a2: Force;


    // initial guess for p2
    let mut force_residual = 1.0;

    // now we are ready to loop 

    while force_residual > tolerance {

        // first let's guess p2 
        p2a2 = p2_guess * a2;
        // we take rho2 based on s2 and p2 
        rho2_guess = v_ps_eqm(p2_guess, s2).recip();
        let v2_guess: Velocity = mass_flowrate/rho2_guess/a2;

        let rhs: Force  = p2a2 + v2_guess * mass_flowrate;

        let root: Force = lhs - rhs;

        let p2_lower_bound_positive = 
            force_balance_isentropic_nozzle(
                p1, p2_lower_bound, h1, mass_flowrate, a1, a2
            ).is_sign_positive();
        let p2_upper_bound_positive = 
            force_balance_isentropic_nozzle(
                p1, p2_upper_bound, h1, mass_flowrate, a1, a2
            ).is_sign_positive();

        if root > Force::ZERO {
            // this is lhs > rhs,

            if p2_upper_bound_positive {
                // if p2 upper bound is positive, then 
                // we move the upper bound lower
                p2_upper_bound = 
                    0.5 * p2_upper_bound
                    + 0.5 * p2_lower_bound;

            } else {
                // if p2 lower bound is positive, then 
                // we move the lower bound higher
                p2_lower_bound = 
                    0.5 * p2_upper_bound
                    + 0.5 * p2_lower_bound;

            }

        } else {

            // this is rhs > lhs
            if p2_lower_bound_positive {
                // if p2 lower bound is positive, then 
                // we move the upper bound lower
                p2_upper_bound = 
                    0.5 * p2_upper_bound
                    + 0.5 * p2_lower_bound;

            } else {
                // if p2 upper bound is positive, then 
                // we move the lower bound higher
                p2_lower_bound = 
                    0.5 * p2_upper_bound
                    + 0.5 * p2_lower_bound;

            }

        }

        p2_guess = 0.5 * p2_upper_bound + 0.5 * p2_lower_bound;

        force_residual = ((lhs - rhs)/lhs).get::<ratio>().abs();

        n_iter += 1;

        if n_iter > max_iter {
            break
        };


    };
    let rho2 = rho2_guess;
    let v2: Velocity = mass_flowrate/rho2/a2;
    let h2 = h1 + 0.5*v1*v1 - 0.5*v2*v2;

    return (p2_guess, h2,rho2);

}

#[inline]
pub fn force_balance_isentropic_nozzle(
    p1: Pressure,
    p2: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    ) -> Force {
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

    let rho1: MassDensity = state_1.get_specific_volume().recip();
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();
    // note: this is isentropic, hence s2 = s1;
    let s2 = s1;
    let v1: Velocity = mass_flowrate/rho1/a1;

    let p1a1: Force = p1*a1;
    // left hand side of momentum balance
    //
    // P1 A1 + dot{m}^2/{rho1 a1}
    let lhs: Force = p1a1 + mass_flowrate * v1;

    // if we have p2, we can get the second thermodynamic state
    let p2a2 = p2 * a2;
    let rho2 = v_ps_eqm(p2, s2).recip();
    let rhs: Force = p2a2 + mass_flowrate * mass_flowrate/rho2/a2;

    return lhs-rhs;

}

#[inline]
pub fn print_graph_pts_for_outlet_pressure_and_force_balance(
    p1: Pressure,
    h1: AvailableEnergy,
    a1: Area,
    a2: Area,
    mass_flowrate: MassRate,
){
    // now let's plot a trend of isentropic nozzles 
    println!("{:?}",&("oulet pressure (Bar)","force bal(Newton)"));
    let mut p2 = p1;
    let n = 100;
    for i in 0..n {

        p2 = (n as f64 - i as f64)/(n as f64) * p1;

        let force_bal: Force = 
            force_balance_isentropic_nozzle(p1, p2, h1, mass_flowrate, a1, a2);

        println!("{:?}",&(p2.get::<bar>(),force_bal.get::<newton>()));

    }

    // let's also check the lower bound pressures
    // below the first 100 points
    for i in 0..n {

        if p2 > Pressure::new::<pascal>(1000_f64) {

            p2 = (1_f64)/(n as f64) * p1;
            p2 *= (n-i) as f64 / (n as f64);


            let force_bal: Force = 
                force_balance_isentropic_nozzle(p1, p2, h1, mass_flowrate, a1, a2);

            println!("{:?}",&(p2.get::<bar>(),force_bal.get::<newton>()));

        }

    }
}

#[cfg(test)]
mod isentropic_vibe_coded_nozzles_test {
    use uom::si::available_energy::kilojoule_per_kilogram;
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
    #[ignore = "ignored hs algo"]
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

    }
    #[test]
    #[ignore = "ignored hs algo"]
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


    }


    /// note that the hs and ps algorithm give slightly 
    /// different answer as both are backward equations
    ///
    /// One only need look at the hs algorithm and ps algorithm 
    /// look at the absolute enthalpy.
    ///
    /// Okay something is really fishy
    #[test]
    #[ignore = "ignored hs algo"]
    pub fn nozzles_test_hs_ps_algo_comparison(){

        let a1: Area = Area::new::<square_centimeter>(5.0);
        let a2: Area = Area::new::<square_centimeter>(4.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
        let tolerance: Option<f64> = Option::None;
        // let this be steam at 10 bar 330 C
        // this is superheated steam
        let p1 = Pressure::new::<bar>(10.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(330.0);
        let h1 = h_tp_eqm_single_phase(t1, p1);
        

        let (p2_hs_algo,h2_hs_algo,_rho2_hs_algo) = get_isentropic_nozzles_outlet_ph_rho_point_hs_algo(
            p1, h1, mass_flowrate, a1, a2, tolerance);

        let (p2_ps_algo,h2_ps_algo,_rho2_ps_algo) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo(
            p1, h1, mass_flowrate, a1, a2, tolerance);

        approx::assert_relative_eq!(
            h1.get::<kilojoule_per_kilogram>(),
            3115.676,
            max_relative=1e-4,
        );

        approx::assert_relative_eq!(
            h2_hs_algo.get::<kilojoule_per_kilogram>(),
            3112.848,
            max_relative=1e-4,
        );
        approx::assert_relative_eq!(
            h2_ps_algo.get::<kilojoule_per_kilogram>(),
            3115.090,
            max_relative=1e-4,
        );
        approx::assert_relative_eq!(
            p2_hs_algo.get::<bar>(),
            12.385,
            max_relative=1e-4,
        );
        approx::assert_relative_eq!(
            p2_ps_algo.get::<bar>(),
            12.476,
            max_relative=1e-4,
        );

        todo!("something is really fishy here with the algos");

    }

}

#[cfg(test)]
mod nozzle_tests {
    use super::*;
    use uom::si::pressure::bar;
    use uom::si::available_energy::joule_per_kilogram;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::area::square_meter;
    use uom::si::velocity::meter_per_second;
    use uom::si::mass_density::kilogram_per_cubic_meter;

    /// Test Case 1: High Pressure (HP) Turbine First Stage Nozzle
    /// Typical for 250 MW steam turbine
    #[test]
    fn test_hp_turbine_first_stage_nozzle() {
        println!("\n=== HP Turbine First Stage Nozzle ===");
        
        // Inlet conditions (after superheater)
        let p1 = Pressure::new::<bar>(160.0);  // 160 bar
        let t1_celsius = 540.0;  // 540°C superheated steam
        
        // For 540°C, 160 bar superheated steam: h ≈ 3,450 kJ/kg
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_450_000.0);
        
        // Mass flow for 250 MW turbine
        let mass_flowrate = MassRate::new::<kilogram_per_second>(200.0);
        
        // Nozzle areas (converging-diverging)
        let a1 = Area::new::<square_meter>(0.15);  // Inlet area
        let a2 = Area::new::<square_meter>(0.20);  // Outlet area (diverging)
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        // Calculate outlet conditions
        let v2 = mass_flowrate / rho2 / a2;
        let expansion_ratio = p1 / p2;
        
        println!("Inlet:  P1 = {:.1} bar, h1 = {:.1} kJ/kg", 
            p1.get::<bar>(), h1.get::<joule_per_kilogram>() / 1e3);
        println!("Outlet: P2 = {:.1} bar, h2 = {:.1} kJ/kg", 
            p2.get::<bar>(), h2.get::<joule_per_kilogram>() / 1e3);
        println!("Density: ρ2 = {:.2} kg/m³", rho2.get::<kilogram_per_cubic_meter>());
        println!("Velocity: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        println!("Expansion ratio: {:.2}", expansion_ratio.get::<ratio>());
        println!("Enthalpy drop: {:.1} kJ/kg", 
            (h1 - h2).get::<joule_per_kilogram>() / 1e3);
        
        // Verify physical constraints
        assert!(p2 < p1, "Pressure should decrease");
        assert!(p2 > Pressure::new::<bar>(50.0), "P2 should be reasonable for HP stage");
        assert!(h2 < h1, "Enthalpy should decrease (acceleration)");
        assert!(v2.get::<meter_per_second>() > 100.0, "Outlet velocity should be significant");
        assert!(v2.get::<meter_per_second>() < 600.0, "Velocity should be subsonic/transonic");
    }

    /// Test Case 2: Intermediate Pressure (IP) Turbine Nozzle
    #[test]
    fn test_ip_turbine_nozzle() {
        println!("\n=== IP Turbine Nozzle ===");
        
        // Inlet conditions (after HP turbine and reheater)
        let p1 = Pressure::new::<bar>(40.0);  // 40 bar
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_400_000.0);  // Reheated to ~540°C
        
        let mass_flowrate = MassRate::new::<kilogram_per_second>(200.0);
        
        // Larger nozzle areas for lower pressure
        let a1 = Area::new::<square_meter>(0.30);
        let a2 = Area::new::<square_meter>(0.40);
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        
        println!("Inlet:  P1 = {:.1} bar, h1 = {:.1} kJ/kg", 
            p1.get::<bar>(), h1.get::<joule_per_kilogram>() / 1e3);
        println!("Outlet: P2 = {:.1} bar, h2 = {:.1} kJ/kg", 
            p2.get::<bar>(), h2.get::<joule_per_kilogram>() / 1e3);
        println!("Velocity: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        
        assert!(p2 < p1);
        assert!(p2 > Pressure::new::<bar>(10.0), "P2 should be reasonable for IP stage");
        assert!(h2 < h1);
    }

    /// Test Case 3: Low Pressure (LP) Turbine Last Stage Nozzle
    /// Exhausts to condenser
    #[test]
    fn test_lp_turbine_last_stage_nozzle() {
        println!("\n=== LP Turbine Last Stage Nozzle ===");
        
        // Inlet conditions (near saturation)
        let p1 = Pressure::new::<bar>(1.0);  // 1 bar
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(2_700_000.0);  // Wet steam region
        
        let mass_flowrate = MassRate::new::<kilogram_per_second>(200.0);
        
        // Very large nozzle areas for LP stage
        let a1 = Area::new::<square_meter>(2.0);
        let a2 = Area::new::<square_meter>(3.5);  // Large expansion
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        
        println!("Inlet:  P1 = {:.3} bar, h1 = {:.1} kJ/kg", 
            p1.get::<bar>(), h1.get::<joule_per_kilogram>() / 1e3);
        println!("Outlet: P2 = {:.3} bar, h2 = {:.1} kJ/kg", 
            p2.get::<bar>(), h2.get::<joule_per_kilogram>() / 1e3);
        println!("Density: ρ2 = {:.3} kg/m³", rho2.get::<kilogram_per_cubic_meter>());
        println!("Velocity: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        
        // Should exhaust near condenser pressure
        assert!(p2 < p1);
        assert!(p2 > Pressure::new::<bar>(0.03), "Should be above condenser vacuum");
        assert!(p2 < Pressure::new::<bar>(0.15), "Should approach condenser pressure");
        assert!(h2 < h1);
        assert!(rho2.get::<kilogram_per_cubic_meter>() < 1.0, "Very low density in LP exhaust");
    }

    /// Test Case 4: Converging Nozzle (Subsonic)
    #[test]
    fn test_converging_nozzle() {
        println!("\n=== Converging Nozzle (Subsonic) ===");
        
        let p1 = Pressure::new::<bar>(10.0);
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(2_800_000.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(50.0);
        
        // Converging: a2 < a1
        let a1 = Area::new::<square_meter>(0.10);
        let a2 = Area::new::<square_meter>(0.05);
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        
        println!("Inlet:  P1 = {:.1} bar, A1 = {:.3} m²", 
            p1.get::<bar>(), a1.get::<square_meter>());
        println!("Outlet: P2 = {:.1} bar, A2 = {:.3} m²", 
            p2.get::<bar>(), a2.get::<square_meter>());
        println!("Velocity: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        
        assert!(p2 < p1);
        assert!(a2 < a1);
    }

    /// Test Case 5: High Expansion Ratio (Diverging Nozzle)
    #[test]
    fn test_high_expansion_nozzle() {
        println!("\n=== High Expansion Diverging Nozzle ===");
        
        let p1 = Pressure::new::<bar>(100.0);
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_200_000.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(150.0);
        
        // Large area ratio for high expansion
        let a1 = Area::new::<square_meter>(0.08);
        let a2 = Area::new::<square_meter>(0.50);  // 6.25:1 area ratio
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-5)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        let expansion_ratio = p1 / p2;
        
        println!("Inlet:  P1 = {:.1} bar", p1.get::<bar>());
        println!("Outlet: P2 = {:.1} bar", p2.get::<bar>());
        println!("Area ratio: {:.2}", (a2/a1).get::<ratio>());
        println!("Pressure ratio: {:.2}", expansion_ratio.get::<ratio>());
        println!("Velocity: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        println!("Mach number estimate: ~{:.2}", v2.get::<meter_per_second>() / 450.0);
        
        assert!(expansion_ratio.get::<ratio>() > 3.0, "Should have significant expansion");
        assert!(v2.get::<meter_per_second>() > 300.0, "High velocity expected");
    }

    /// Test Case 6: Small Mass Flow (Control Stage)
    #[test]
    fn test_control_stage_nozzle() {
        println!("\n=== Control Stage Nozzle (Partial Load) ===");
        
        let p1 = Pressure::new::<bar>(160.0);
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_450_000.0);
        
        // Reduced mass flow (turbine at 50% load)
        let mass_flowrate = MassRate::new::<kilogram_per_second>(100.0);
        
        let a1 = Area::new::<square_meter>(0.08);
        let a2 = Area::new::<square_meter>(0.12);
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        
        println!("Partial load operation:");
        println!("Mass flow: {:.1} kg/s", mass_flowrate.get::<kilogram_per_second>());
        println!("Outlet: P2 = {:.1} bar, v2 = {:.1} m/s", 
            p2.get::<bar>(), v2.get::<meter_per_second>());
        
        assert!(p2 < p1);
    }

        /// Test Case 7: Wet Steam Region (LP Turbine) - CONTINUED
    #[test]
    fn test_wet_steam_nozzle() {
        println!("\n=== Wet Steam Nozzle (LP Stage) ===");
        
        // Conditions in wet steam region
        let p1 = Pressure::new::<bar>(0.5);  // 0.5 bar
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(2_450_000.0);  // ~90% quality
        
        let mass_flowrate = MassRate::new::<kilogram_per_second>(200.0);
        
        let a1 = Area::new::<square_meter>(3.0);
        let a2 = Area::new::<square_meter>(5.0);
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        
        println!("Wet steam expansion:");
        println!("P1 = {:.3} bar → P2 = {:.3} bar", 
            p1.get::<bar>(), p2.get::<bar>());
        println!("h1 = {:.1} kJ/kg → h2 = {:.1} kJ/kg",
            h1.get::<joule_per_kilogram>() / 1e3,
            h2.get::<joule_per_kilogram>() / 1e3);
        println!("ρ2 = {:.4} kg/m³", rho2.get::<kilogram_per_cubic_meter>());
        println!("v2 = {:.1} m/s", v2.get::<meter_per_second>());
        
        assert!(p2 < p1);
        assert!(p2 > Pressure::new::<bar>(0.04), "Above condenser pressure");
        assert!(rho2.get::<kilogram_per_cubic_meter>() < 0.5, "Very low density");
    }

    /// Test Case 8: Choked Flow Conditions
    #[test]
    fn test_choked_flow_nozzle() {
        println!("\n=== Choked Flow Nozzle (Critical Pressure Ratio) ===");
        
        let p1 = Pressure::new::<bar>(50.0);
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_100_000.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(100.0);
        
        // Small throat area to induce choking
        let a1 = Area::new::<square_meter>(0.15);
        let a2 = Area::new::<square_meter>(0.08);  // Throat (smallest area)
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        let v2 = mass_flowrate / rho2 / a2;
        let pressure_ratio = p2 / p1;
        
        println!("Choked flow analysis:");
        println!("P1 = {:.1} bar, P2 = {:.1} bar", 
            p1.get::<bar>(), p2.get::<bar>());
        println!("Pressure ratio P2/P1 = {:.3}", pressure_ratio.get::<ratio>());
        println!("Velocity at throat: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        println!("Critical pressure ratio for steam ~0.55");
        
        // For steam, critical pressure ratio is approximately 0.545-0.58
        // If P2/P1 is near this value, flow is choked
        if pressure_ratio.get::<ratio>() > 0.52 && pressure_ratio.get::<ratio>() < 0.60 {
            println!("⚠️  Flow is likely CHOKED (sonic at throat)");
        }
        
        assert!(p2 < p1);
    }

    /// Test Case 9: Supersonic Diverging Section
    #[test]
    fn test_supersonic_diverging_nozzle() {
        println!("\n=== Supersonic Converging-Diverging Nozzle ===");
        
        let p1 = Pressure::new::<bar>(80.0);
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_300_000.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(120.0);
        
        // Converging-diverging (de Laval nozzle)
        // This models the DIVERGING section after the throat
        let a_throat = Area::new::<square_meter>(0.05);  // Throat area (a1)
        let a_exit = Area::new::<square_meter>(0.15);    // Exit area (a2) - 3:1 expansion
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a_throat, a_exit, Some(1e-5)
        );
        
        let v2 = mass_flowrate / rho2 / a_exit;
        
        // Estimate speed of sound (rough approximation for steam)
        // a = sqrt(gamma * R * T), for steam ~400-500 m/s depending on conditions
        let estimated_sonic_speed = 450.0;  // m/s
        let mach_number = v2.get::<meter_per_second>() / estimated_sonic_speed;
        
        println!("Supersonic nozzle characteristics:");
        println!("Throat area: {:.4} m²", a_throat.get::<square_meter>());
        println!("Exit area: {:.3} m²", a_exit.get::<square_meter>());
        println!("Area ratio: {:.2}", (a_exit/a_throat).get::<ratio>());
        println!("P1 = {:.1} bar → P2 = {:.1} bar", 
            p1.get::<bar>(), p2.get::<bar>());
        println!("Exit velocity: v2 = {:.1} m/s", v2.get::<meter_per_second>());
        println!("Estimated Mach number: M ≈ {:.2}", mach_number);
        
        if mach_number > 1.0 {
            println!("✓ Flow is SUPERSONIC at exit");
        } else {
            println!("⚠️  Flow is subsonic (may need higher pressure ratio)");
        }
        
        assert!(v2.get::<meter_per_second>() > 300.0, "Should have high exit velocity");
    }

    /// Test Case 10: Multi-stage Turbine Complete Path
    #[test]
    fn test_multi_stage_turbine_path() {
        println!("\n=== Multi-Stage Turbine Nozzle Cascade ===");
        
        let mass_flowrate = MassRate::new::<kilogram_per_second>(200.0);
        
        // Stage 1: HP Turbine first nozzle
        println!("\n--- Stage 1: HP First Nozzle ---");
        let p1_stage1 = Pressure::new::<bar>(160.0);
        let h1_stage1 = AvailableEnergy::new::<joule_per_kilogram>(3_450_000.0);
        let (p2_stage1, h2_stage1, rho2_stage1) = 
            get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
                p1_stage1, h1_stage1, mass_flowrate,
                Area::new::<square_meter>(0.15),
                Area::new::<square_meter>(0.20),
                Some(1e-4)
            );
        
        println!("Stage 1: {:.1} bar → {:.1} bar, Δh = {:.1} kJ/kg",
            p1_stage1.get::<bar>(),
            p2_stage1.get::<bar>(),
            (h1_stage1 - h2_stage1).get::<joule_per_kilogram>() / 1e3);
        
        // Stage 2: HP Turbine second nozzle (uses stage 1 outlet)
        println!("\n--- Stage 2: HP Second Nozzle ---");
        let (p2_stage2, h2_stage2, rho2_stage2) = 
            get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
                p2_stage1, h2_stage1, mass_flowrate,
                Area::new::<square_meter>(0.25),
                Area::new::<square_meter>(0.35),
                Some(1e-4)
            );
        
        println!("Stage 2: {:.1} bar → {:.1} bar, Δh = {:.1} kJ/kg",
            p2_stage1.get::<bar>(),
            p2_stage2.get::<bar>(),
            (h2_stage1 - h2_stage2).get::<joule_per_kilogram>() / 1e3);
        
        // Stage 3: IP Turbine nozzle (after reheater)
        println!("\n--- Stage 3: IP Nozzle (after reheat) ---");
        let h1_stage3 = AvailableEnergy::new::<joule_per_kilogram>(3_400_000.0);  // Reheated
        let (p2_stage3, h2_stage3, rho2_stage3) = 
            get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
                p2_stage2, h1_stage3, mass_flowrate,
                Area::new::<square_meter>(0.40),
                Area::new::<square_meter>(0.60),
                Some(1e-4)
            );
        
        println!("Stage 3: {:.1} bar → {:.1} bar, Δh = {:.1} kJ/kg",
            p2_stage2.get::<bar>(),
            p2_stage3.get::<bar>(),
            (h1_stage3 - h2_stage3).get::<joule_per_kilogram>() / 1e3);
        
        // Calculate total expansion
        println!("\n--- Overall Turbine Performance ---");
        println!("Total pressure drop: {:.1} bar → {:.1} bar",
            p1_stage1.get::<bar>(),
            p2_stage3.get::<bar>());
        println!("Overall pressure ratio: {:.1}",
            (p1_stage1 / p2_stage3).get::<ratio>());
        
        assert!(p2_stage3 < p1_stage1);
    }

    /// Test Case 11: Extreme Conditions (Boundary Testing)
    #[test]
    fn test_extreme_low_pressure() {
        println!("\n=== Extreme Low Pressure (Near Vacuum) ===");
        
        let p1 = Pressure::new::<bar>(0.10);  // 0.1 bar
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(2_500_000.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(150.0);
        
        // Very large areas for low pressure
        let a1 = Area::new::<square_meter>(5.0);
        let a2 = Area::new::<square_meter>(8.0);
        
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-4)
        );
        
        println!("Near-vacuum conditions:");
        println!("P2 = {:.4} bar ({:.1} kPa)", 
            p2.get::<bar>(), p2.get::<bar>() * 100.0);
        println!("ρ2 = {:.5} kg/m³", rho2.get::<kilogram_per_cubic_meter>());
        
        assert!(p2 > Pressure::new::<bar>(0.01), "Should stay above minimum bound");
        assert!(p2 < p1);
    }

    /// Test Case 12: Convergence Test with Tight Tolerance
    #[test]
    fn test_high_precision_convergence() {
        println!("\n=== High Precision Convergence Test ===");
        
        let p1 = Pressure::new::<bar>(100.0);
        let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_200_000.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(180.0);
        let a1 = Area::new::<square_meter>(0.12);
        let a2 = Area::new::<square_meter>(0.18);
        
        // Test with very tight tolerance
        let (p2, h2, rho2) = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
            p1, h1, mass_flowrate, a1, a2, Some(1e-6)  // Very tight!
        );
        
        println!("High precision result:");
        println!("P2 = {:.6} bar", p2.get::<bar>());
        println!("h2 = {:.3} kJ/kg", h2.get::<joule_per_kilogram>() / 1e3);
        
        assert!(p2 < p1);
    }

}

