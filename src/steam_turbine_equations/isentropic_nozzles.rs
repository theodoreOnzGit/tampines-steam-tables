use uom::ConstZero;
use uom::si::f64::*;
use uom::si::pressure::{bar, millibar};
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

macro_rules! function_debug {
    () => {{
        fn f() {}
        fn type_name_of<T>(_: T) -> &'static str {
            std::any::type_name::<T>()
        }
        let name = type_name_of(f);
        // Gets the function name from the full path
        name.strip_suffix("::f")
            .unwrap_or(name)
            .rsplit("::")
            .next()
            .unwrap_or(name)
    }};
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
    let mut rho2_guess: MassDensity = rho1 * 0.9;
    let mut p2_upper_bound = p1;
    // don't expect turbine exhaust to be 
    // lower than condenser pressure on a good day, 0.04 bar
    let mut p2_lower_bound = Pressure::new::<bar>(0.01);

    // with these two bounds, let's calculate the signs

    fn lhs_minus_rhs_greater_than_zero(
        lhs: Force,
        p2: Pressure,
        s2: SpecificHeatCapacity,
        mass_flowrate: MassRate,
        a2: Area,) -> bool {

        let p2a2 = p2 * a2;
        let rho2 = v_ps_eqm(p2, s2).recip();
        let rhs: Force = p2a2 + mass_flowrate * mass_flowrate/rho2/a2;

        if lhs > rhs {
            return true;
        } else {
            return false;
        }
    }

    let p2_upper_bound_positive = 
        lhs_minus_rhs_greater_than_zero(
            lhs, p2_upper_bound, s2, mass_flowrate, a2
        );
    let p2_lower_bound_positive = 
        lhs_minus_rhs_greater_than_zero(
            lhs, p2_lower_bound, s2, mass_flowrate, a2
        );

    if p2_lower_bound_positive && p2_upper_bound_positive ||
    !p2_lower_bound_positive && !p2_upper_bound_positive {
        dbg!("Error in {}: ", function_debug!());
        panic!("no sign chg, bisection algorithm won't work");
    };



    let mut p2a2: Force;


    // initial guess for p2
    let mut p2_guess: Pressure = 0.5*p2_upper_bound;
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
    let h2 = h1 - 0.5* v2*v2;

    return (p2_guess, h2,rho2);

}



#[cfg(test)]
mod nozzles_test {
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
