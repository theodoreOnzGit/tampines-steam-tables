use uom::ConstZero;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::pressure::megapascal;
use uom::si::pressure::pascal;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::volume::cubic_meter;

use crate::constants::p_crit_water;
use crate::constants::t_crit_water;
use crate::interfaces::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::interfaces::functional_programming::ph_flash_eqm::s_ph_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::t_ph_eqm;
use crate::interfaces::functional_programming::ps_flash_eqm::ps_flash_region;
use crate::interfaces::functional_programming::ps_flash_eqm::w_ps_wood_wallis;
use crate::interfaces::functional_programming::ps_flash_eqm::x_ps_flash;
use crate::prelude::functional_programming::ps_flash_eqm::mass_flux_ps_eqm_throat;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::w_ph_wood_wallis;
use crate::prelude::functional_programming::ph_flash_eqm::cv_ph_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::cp_ph_eqm;
use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::prelude::functional_programming::pt_flash_eqm::FwdEqnRegion;

use crate::prelude::TampinesSteamTableCV;
use crate::region_2_vapour::s_tp_2;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;


/// This is an algorithm to obtain outlet thermodynamic state 
/// for a converging nozzle with subsonic flow
/// Given inlet conditions (p1, h1, v1) and geometry (a1, a2), calculates
/// the throat conditions assuming choked flow (M = 1 at exit).
///
/// # Arguments
/// * `p1` - Inlet pressure
/// * `h1` - Inlet specific enthalpy
/// * `a1` - Inlet area
/// * `a2` - Throat (exit) area
/// * `v1` - Inlet velocity
///
/// # Returns
/// Tuple of (p2, h2, mass_flowrate) where:
/// * `p2` - Throat pressure
/// * `h2` - Throat specific enthalpy
/// * `mass_flowrate` - Choked mass flow rate (determined by throat conditions)
///
/// # Warnings
/// - Warns if mass balance error > 5% (inlet vs throat)
/// - Warns if momentum balance error > 5%
#[inline]
pub fn get_choked_flow_state_for_nozzle_subsonic_to_sonic(
    p1: Pressure,
    h1: AvailableEnergy,
    a1: Area,
    a2: Area,
    v1: Velocity,
) -> (Pressure, AvailableEnergy, MassRate) {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let rho_1 = state_1.get_rho();


    let s1 = state_1.get_specific_entropy();

    let h0 = h1 + 0.5 * v1 * v1;
    // stagnation state (initial)
    let state_0 = TampinesSteamTableCV::new_from_hs(h0, s1, ref_vol);
    let p0 = state_0.get_pressure();


    // now, i'll have to get a solver for choked flow 

    // let's use the critical pressure 

    let critical_pressure_ratio: Ratio = 
        state_0.get_critical_pressure_ratio_pure_vapour();

    // this is critical pressure for mach 1
    let p2 = critical_pressure_ratio * p0;
    // let's get speed of sound here 
    let s2 = s1;
    let state_2 = TampinesSteamTableCV::new_from_ps(p2, s2, ref_vol);
    let c = state_2.get_speed_of_sound();
    let h2 = state_2.get_specific_enthalpy();
    let rho_2 = state_2.get_rho();

    // now let's get mass balance first 

    let mass_flowrate = a1 * v1 * rho_1;
    let choked_mass_flowrate = a2 * c * rho_2;

    // assert that it isn't too different 
    //

    let mass_flowrate_error: Ratio = 
        (choked_mass_flowrate - mass_flowrate)/mass_flowrate;

    if mass_flowrate_error.abs() >= Ratio::new::<ratio>(0.05) {
        eprintln!("Warning: Mass balance error = {:.2}%. \
                   Inlet conditions may not be consistent with choked flow.",
                   mass_flowrate_error.get::<ratio>() * 100.0);
    }
    
    let momentum_balance_check_1d = false; 
    
    // the mommentum balance check here does not account for wall forces,
    // so it may not always work
    //
    // so for isentropic flows, we only do mass and energy balances 
    // the force balances are not accounted for
    if momentum_balance_check_1d {

        // after this, it is time to do force balance
        // force balance is 
        //
        // P1A1 + v1 * mass_flowrate  = P2A2 + c*mass_flowrate

        let p1a1: Force = p1 * a1;
        let p2a2: Force = p2 * a2;

        let lhs: Force = p1a1 + v1 * choked_mass_flowrate;
        let rhs: Force = p2a2 + c * choked_mass_flowrate;

        // for nozzle, 
        // mass flowrate =  (p1a1 - p2a2) / (-v1 + c)

        let momentum_balance_error: Ratio = (lhs - rhs)/lhs;

        if momentum_balance_error.abs() >= Ratio::new::<ratio>(0.05) {
            eprintln!("Warning: Momentum balance error = {:.2}%. \
                   Check if flow is truly choked.",
                   momentum_balance_error.get::<ratio>() * 100.0);
            // note: example question didn't do this
            //
            // but to be fair, the nozzle walls exert a force on the fluid, 
            // this is note accounted for in the above momentum balance
            dbg!(&(lhs,rhs));
        }
    }

    
    return (p2, h2, choked_mass_flowrate);

}

/// based on steam pressure and enthalpy, 
/// as well as mass flowrate, obtain the choked flow area
/// you have to give the stagnation enthalpy and pressure as 
/// inputs
pub fn get_choked_flow_nozzle_area(
    p0: Pressure,
    h0: AvailableEnergy,
    mass_flowrate: MassRate,
    ) -> Area {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_0 = TampinesSteamTableCV::new_from_ph(p0, h0, ref_vol);
    let s0 = state_0.get_specific_entropy();
    let s1 = s0;


    // now, i'll have to get a solver for choked flow 

    // let's use the critical pressure 


    let critical_pressure_ratio: Ratio = 
        state_0.get_critical_pressure_ratio_pure_vapour();

    // this is critical pressure for mach 1
    let p2 = critical_pressure_ratio * p0;
    // let's get speed of sound here 
    let s2 = s1;
    let state_2 = TampinesSteamTableCV::new_from_ps(p2, s2, ref_vol);
    let c = state_2.get_speed_of_sound();
    let rho_2 = state_2.get_rho();

    // after getting c, we should be able to get area
    let a_throat: Area = mass_flowrate/rho_2/c;

    return a_throat;
    
}
/// based on steam pressure and enthalpy, 
/// as well as mass flowrate, obtain the choked flow area
/// you have to give the stagnation enthalpy and pressure as 
/// inputs
///
/// isentropic nozzle assumed
///
/// also, mass flowrate and exit pressure given,
/// assuming throat velocity is speed of sound
pub fn get_choked_flow_supersonic_nozzle_exit_area_and_state(
    p_throat: Pressure,
    h_throat: AvailableEnergy,
    p_exit: Pressure,
    mass_flowrate: MassRate,
    ) -> (Pressure, AvailableEnergy, Area) {


    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_throat = TampinesSteamTableCV::new_from_ph(
        p_throat, h_throat, ref_vol
    );

    let c = state_throat.get_speed_of_sound();
    // stagnation properties
    let s0 = state_throat.get_specific_entropy();
    // I assume speed of sound at the throat
    let h0 = h_throat + 0.5 * c * c;

    // at exit, the pressure is given

    let state_exit = TampinesSteamTableCV::new_from_ps(
        p_exit, s0, ref_vol
    );

    let h_exit = state_exit.get_specific_enthalpy();
    let v_exit: Velocity = (2.0 * (h0 - h_exit)).sqrt();
    let rho_exit = state_exit.get_rho();
    let a_exit = mass_flowrate/v_exit/rho_exit;

    


    return (p_exit, h_exit, a_exit);
    
}

#[cfg(test)]
mod choked_flow_examples{
    use uom::si::area::square_centimeter;
    use uom::si::available_energy::kilojoule_per_kilogram;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::{kilopascal, megapascal};
    use uom::si::f64::*;
    use uom::si::ratio::ratio;
    use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::velocity::meter_per_second;
    use uom::si::volume::cubic_meter;

    use crate::prelude::TampinesSteamTableCV;
    use crate::steam_turbine_equations::choked_flow::{get_choked_flow_supersonic_nozzle_exit_area_and_state, get_choked_flow_nozzle_area, get_choked_flow_state_for_nozzle_subsonic_to_sonic};


    /// this is from example 17-16 in Cengel's thermodynamics 8th edition
    ///
    /// Steam enters a converging-diverging nozzle at 2 MPa and 400C 
    /// with negligble velocity and flowrate is 2.5 kg/s 
    /// it exits at pressure of 300 kPa 
    ///
    /// from nozzle entrance to throat, it is isentropic 
    /// after throat, it is 93 % efficient 
    ///
    /// what is 
    /// 1. the throat and exit area
    /// 2. mach number at throat and nozzle exit
    #[test]
    pub fn steam_flow_cd_nozzle(){

        let p1 = Pressure::new::<megapascal>(2.0);
        let t1 = ThermodynamicTemperature::new::<degree_celsius>(400.0);
        let mass_flowrate = MassRate::new::<kilogram_per_second>(2.5);
        let p_exit = Pressure::new::<kilopascal>(300.0);

        // stagnation pressure = p1 as velocity is negligble 

        let p0 = p1;
        // we assume steam is superheated, so quality is 1 
        let ref_vol = Volume::new::<cubic_meter>(1.0);

        let state_0 = TampinesSteamTableCV::new_from_tp_quality_1(
            t1, p0, ref_vol
        );

        let h0 = state_0.get_specific_enthalpy();

        // let's get entropy
        let s0 = state_0.get_specific_entropy();

        let critical_pressure_ratio = 
            state_0.get_critical_pressure_ratio_pure_vapour();

        // In Cengel's textbook, critical pressure ratio is approximated 
        // as 0.546 

        approx::assert_relative_eq!(
            critical_pressure_ratio.get::<ratio>(),
            0.546,
            max_relative=1e-3,
        );

        // critical pressure ratio checks out!
        // now, we expect throat pressure to be 1.09 MPa
        let p_throat = p1 * critical_pressure_ratio;
        approx::assert_relative_eq!(
            p_throat.get::<megapascal>(),
            1.09,
            max_relative=5e-3,
        );

        // now for throat conditions 

        let state_throat = 
            TampinesSteamTableCV::new_from_ps(
                p_throat, s0, ref_vol
            );

        // let's find the velocity here 
        // this is, of course, the speed of sound
        let c = state_throat.get_speed_of_sound();

        // now, based on enthalpy balance, we should get 
        // a h_throat 

        let h_throat = state_throat.get_specific_enthalpy();

        approx::assert_relative_eq!(
            h_throat.get::<kilojoule_per_kilogram>(),
            3076.8,
            max_relative=1e-3,
        );

        let v_throat: Velocity = (2.0 *(h0-h_throat)).sqrt();

        // if we can see here, v_throat and c is about the same
        // Cengel's answer is 585.8 m/s
        // which is about there
        approx::assert_relative_eq!(
            c.get::<meter_per_second>(),
            584.74,
            max_relative=1e-3,
        );

        approx::assert_relative_eq!(
            v_throat.get::<meter_per_second>(),
            584.74,
            max_relative=1e-3,
        );

        // now, we find throat area is 10.33 cm^2 from Cengel's 
        // that's the answer, but we will work backwards using our function

        // velocity at the inlet is based on the mass flowrate 

        let a2 = Area::new::<square_centimeter>(10.33);
        // now, this a1 is arbitrary.., we only need know that v1 is negligble
        let a1 = 100.0 * a2;
        let rho_1 = state_0.get_rho();
        let v1: Velocity = mass_flowrate/rho_1/a1;
        
        approx::assert_relative_eq!(
            v1.get::<meter_per_second>(),
            3.659,
            max_relative=1e-3,
        );

        // let's see whether given these areas, we get the correct
        // critical pressure predictions 

        let (p_crit, h_crit, mass_flowrate) = 
            get_choked_flow_state_for_nozzle_subsonic_to_sonic(
                p1, h0, a1, a2, v1
            );

        approx::assert_relative_eq!(
            p_crit.get::<megapascal>(),
            1.092,
            max_relative=1e-3,
        );
        approx::assert_relative_eq!(
            h_crit.get::<kilojoule_per_kilogram>(),
            3076.8,
            max_relative=1e-3,
        );
        approx::assert_relative_eq!(
            mass_flowrate.get::<kilogram_per_second>(),
            2.5,
            max_relative=1e-3,
        );

        dbg!(&(
                critical_pressure_ratio,
                p_throat.get::<megapascal>(),
                v_throat.get::<meter_per_second>(),
                c.get::<meter_per_second>(),
        ));

        // now we want to deal with exit properties,
        // we know the pressure and 
        // Cengel writes the entahlpy as 2816.1 kJ/kg given a 0.93 
        // efficiency (we are not testing for that here)

        let h_exit = AvailableEnergy::new::<kilojoule_per_kilogram>(2816.1);

        let state_exit = 
            TampinesSteamTableCV::new_from_ph(
                p_exit, h_exit, ref_vol
            );

        let c_exit: Velocity = state_exit.get_speed_of_sound();


        // from Cengel's the exit speed of sound is 515.4 m/s
        approx::assert_relative_eq!(
            c_exit.get::<meter_per_second>(),
            515.4,
            max_relative=1e-3,
        );

        // from Cengel's the exit entropy is 7.2019 Kj/kg k

        let s_exit = state_exit.get_specific_entropy();

        approx::assert_relative_eq!(
            s_exit.get::<kilojoule_per_kilogram_kelvin>(),
            7.2019,
            max_relative=1e-3,
        );

        // for exit velocity, it is convenient to use a stagnation enthalpy 
        // as the reference, since it is the same through the 
        // whole nozzle 

        let v_exit = (2.0*(h0 - h_exit)).sqrt();
        approx::assert_relative_eq!(
            v_exit.get::<meter_per_second>(),
            929.8,
            max_relative=1e-3,
        );

        // exit mach number is 1.804 

        let mach_number_exit = state_exit.get_mach_number(v_exit);
        approx::assert_relative_eq!(
            mach_number_exit.get::<ratio>(),
            1.804,
            max_relative=1e-3,
        );

        dbg!(&(
                p_exit.get::<megapascal>(),
                v_throat.get::<meter_per_second>(),
                c_exit.get::<meter_per_second>(),
                v_exit.get::<meter_per_second>(),
                s_exit.get::<kilojoule_per_kilogram_kelvin>(),
                mach_number_exit.get::<ratio>(),
        ));


    }

    /// this is a converging diverging nozzle using Cengel's question 
    /// 17-110
    #[test]
    fn steam_flow_cd_nozzle_2(){

        let inlet_temperature = 
            ThermodynamicTemperature::new::<degree_celsius>(
                500.0
            );

        let inlet_pressure = 
            Pressure::new::<megapascal>(
                1.0 
            );

        let mass_flowrate = MassRate::new::<kilogram_per_second>(
            2.5
        );

        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_0 = TampinesSteamTableCV::new_from_tp_quality_1(
            inlet_temperature, inlet_pressure, ref_vol
        );
        let h0 = state_0.get_specific_enthalpy();
        let p0 = inlet_pressure;

        // note, throat diameter is 22.4 cm^2
        let area = get_choked_flow_nozzle_area(p0, h0, mass_flowrate);
        approx::assert_relative_eq!(
            area.get::<square_centimeter>(),
            22.399,
            max_relative=1e-3,
        );
        let ref_vol = Volume::new::<cubic_meter>(1.0);
        let state_0 = TampinesSteamTableCV::new_from_ph(p0, h0, ref_vol);
        let s0 = state_0.get_specific_entropy();
        let s1 = s0;


        // now, i'll have to get a solver for choked flow 

        // let's use the critical pressure 


        let critical_pressure_ratio: Ratio = 
            state_0.get_critical_pressure_ratio_pure_vapour();

        // this is critical pressure for mach 1
        let p_throat = critical_pressure_ratio * p0;
        // let's get speed of sound here 
        let s2 = s1;
        let state_throat = TampinesSteamTableCV::new_from_ps(p_throat, s2, ref_vol);
        let h_throat = state_throat.get_specific_enthalpy();
        let p_exit = Pressure::new::<kilopascal>(200.0);

        let (p_exit, h_exit, exit_area) = 
            get_choked_flow_supersonic_nozzle_exit_area_and_state(
            p_throat, h_throat, p_exit, mass_flowrate
        );

        approx::assert_relative_eq!(
            exit_area.get::<square_centimeter>(),
            31.5,
            max_relative=1e-3,
        );

        // lastly, get mach number for exit
        
        let v_exit = (2.0*(h0-h_exit)).sqrt();
        let state_exit = TampinesSteamTableCV::new_from_ph(
            p_exit, h_exit, ref_vol
        );
        let mach_number = state_exit.get_mach_number(v_exit);


        // cengel reports 1.738
        // this is within 1% of the value given
        approx::assert_relative_eq!(
            mach_number.get::<ratio>(),
            1.738,
            max_relative=1e-2,
        );
        approx::assert_relative_eq!(
            mach_number.get::<ratio>(),
            1.728,
            max_relative=1e-3,
        );
    }

}


/// gets critical pressure and mass flux for water and steam 
/// given stagnation properties,
/// should work for all given regions of steam table
///
/// Note, this relies on the ph, ps algorithm, for region 5, it is 
/// not thoroughly implemented yet
/// gets critical pressure and mass flux for water and steam 
/// given stagnation properties,
/// should work for all given regions of steam table
///
/// Note, this relies on the ph, ps algorithm, for region 5, it is 
/// not thoroughly implemented yet
/// gets critical pressure and mass flux for water and steam 
/// given stagnation properties,
/// should work for all given regions of steam table
///
/// Note, this relies on the ph, ps algorithm, for region 5, it is 
/// not thoroughly implemented yet
#[inline]
pub fn get_critical_pressure_and_mass_flux_with_stagnation_props(
    s0: SpecificHeatCapacity,
    h0: AvailableEnergy,
    p0: Pressure) -> (Pressure, MassFlux) {

    let debug = true;

    // now before anything, we want to get the region of the scans 
    let region_stagnation_props = ph_flash_region(p0, h0);

    // for high entropy states (s0 >= 9.2 kJ/kg/K), 
    // isentropic depressurisation stays in single phase vapour
    // so we can use the pure vapour algorithm directly
    if s0 >= SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(9.2) {
        if debug { println!("using vapour only algorithm, s0 >= 9.2 kJ/kg/K"); }

        let s0_opt = Some(s0);
        let critical_pressure_choked_flow 
            = get_critical_pressure_pure_vapour_ph_stagnation_properties(
                p0, h0, s0_opt
            );

        let c = w_ps_wood_wallis(critical_pressure_choked_flow, s0);
        let rho_throat = v_ps_eqm(critical_pressure_choked_flow, s0).recip();
        let critical_mass_flux = c * rho_throat;

        if debug {
            dbg!(&(s0, critical_pressure_choked_flow, critical_mass_flux));
        }

        return (critical_pressure_choked_flow, critical_mass_flux);
    }

    // for region 2 superheated vapour with s0 < 9.2 kJ/kg/K,
    // check if isentrope crosses into Region 4 during depressurisation
    // by comparing s0 with s_g at stagnation pressure
    //
    // if s0 >= s_g(p0), isentrope stays in vapour — use pure vapour algorithm
    // if s0 < s_g(p0), isentrope crosses into Region 4 — use generalised algorithm
    if region_stagnation_props == FwdEqnRegion::Region2 {
        let p_crit_water_val = p_crit_water();
        if p0 < p_crit_water_val {
            let t_sat = sat_temp_4(p0);
            let s_g = s_tp_2(t_sat, p0);
            if s0 >= s_g {
                if debug { println!("Region 2, s0 >= s_g, using pure vapour algorithm"); }
                let s0_opt = Some(s0);
                let critical_pressure_choked_flow 
                    = get_critical_pressure_pure_vapour_ph_stagnation_properties(
                        p0, h0, s0_opt
                    );
                let c = w_ps_wood_wallis(critical_pressure_choked_flow, s0);
                let rho_throat = v_ps_eqm(critical_pressure_choked_flow, s0).recip();
                let critical_mass_flux = c * rho_throat;
                if debug {
                    dbg!(&(s0, s_g, critical_pressure_choked_flow, critical_mass_flux));
                }
                return (critical_pressure_choked_flow, critical_mass_flux);
            }
            if debug { println!("Region 2, s0 < s_g, using generalised algorithm"); }
        }
    }

    // for region 3, use vapour algorithm only if above critical temperature
    if region_stagnation_props == FwdEqnRegion::Region3 {
        let t_crit = t_crit_water();
        let t0 = t_ph_eqm(p0, h0);
        if t0 > t_crit {
            if debug { println!("Region 3, T > T_crit, using pure vapour algorithm"); }
            let s0_opt = Some(s0);
            let critical_pressure_choked_flow 
                = get_critical_pressure_pure_vapour_ph_stagnation_properties(
                    p0, h0, s0_opt
                );
            let c = w_ps_wood_wallis(critical_pressure_choked_flow, s0);
            let rho_throat = v_ps_eqm(critical_pressure_choked_flow, s0).recip();
            let critical_mass_flux = c * rho_throat;
            if debug {
                dbg!(&(s0, critical_pressure_choked_flow, critical_mass_flux));
            }
            return (critical_pressure_choked_flow, critical_mass_flux);
        }
        if debug { println!("Region 3, T <= T_crit, using generalised algorithm"); }
    }

    // generalised algorithm for:
    // - Region 1 (subcooled liquid)
    // - Region 4 (wet steam)
    // - Region 2 with s0 < s_g (will flash during depressurisation)
    // - Region 3 liquid-like (T <= T_crit)
    //
    // root finding: G_energy = G_HEM
    // f(p) = G_energy(p) - G_HEM(p) = 0
    // positive when g_energy > g_hem (past choked point, too low pressure)
    // negative when g_hem > g_energy (not yet choked, too high pressure)

    let p_min_steam_table = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    let p_decrement: Pressure = 0.05 * (p0 - p_min_steam_table);

    let root_fn = |p_test_pa: f64| -> f64 {
        let p_test = Pressure::new::<pascal>(p_test_pa);
        let region = ps_flash_region(p_test, s0);

        // no choked flow in subcooled liquid
        if region == FwdEqnRegion::Region1 {
            return -1e10_f64;
        }

        let h_test = h_ps_eqm(p_test, s0);
        let kinetic_energy = h0 - h_test;

        // pressure too low, expansion exceeded stagnation enthalpy
        if kinetic_energy < AvailableEnergy::ZERO {
            return 1e10_f64;
        }

        let rho = v_ps_eqm(p_test, s0).recip();
        let g_energy: MassFlux = rho * (2.0 * kinetic_energy).sqrt();
        let g_hem: MassFlux = mass_flux_ps_eqm_throat(p_test, s0);

        if debug {
            let quality = x_ps_flash(p_test, s0);
            dbg!(&(p_test, region, quality, g_energy, g_hem));
        }

        (g_energy - g_hem).get::<kilogram_per_square_meter_second>()
    };

    // scan downward from p0 to find sign change bracket
    let mut p_scan = p0;
    let mut f_prev = root_fn(p_scan.get::<pascal>());
    let mut p_upper_bracket = p0;
    let mut p_lower_bracket = p_min_steam_table;
    let mut found_bracket = false;

    while (p_scan - p_decrement) > p_min_steam_table {
        p_scan -= p_decrement;
        let f_curr = root_fn(p_scan.get::<pascal>());

        if debug {
            dbg!(&(p_scan, f_prev, f_curr));
        }

        // sign change detected — bracket found
        if f_prev * f_curr < 0.0 {
            p_upper_bracket = p_scan + p_decrement;
            p_lower_bracket = p_scan;
            found_bracket = true;
            if debug {
                println!("bracket found!");
                dbg!(&(p_upper_bracket, p_lower_bracket));
            }
            break;
        }

        f_prev = f_curr;
    }

    if !found_bracket {
        panic!("unable to find bracket for critical mass flux root finding");
    }

    // in-house regula falsi (false position) method
    // converges faster than bisection for smooth functions
    let max_iterations = 100;
    let tolerance = 1e-6_f64; // kg/(m²·s)

    let mut p_low = p_lower_bracket;
    let mut p_high = p_upper_bracket;
    let mut f_low = root_fn(p_low.get::<pascal>());
    let mut f_high = root_fn(p_high.get::<pascal>());
    let mut p_crit = p_low;

    for _ in 0..max_iterations {
        let dp = p_high - p_low;
        let df = f_high - f_low;

        // guard against division by zero
        if df.abs() < 1e-30 {
            break;
        }

        // regula falsi interpolation:
        // p_new = p_low - f_low * (p_high - p_low) / (f_high - f_low)
        p_crit = p_low - Pressure::new::<pascal>(
            f_low * dp.get::<pascal>() / df
        );

        // clamp to bounds
        if p_crit < p_low { p_crit = p_low; }
        if p_crit > p_high { p_crit = p_high; }

        let f_crit = root_fn(p_crit.get::<pascal>());

        if debug {
            dbg!(&(p_crit, f_crit));
        }

        // check convergence
        if f_crit.abs() < tolerance {
            break;
        }

        // update bracket
        if f_low * f_crit < 0.0 {
            p_high = p_crit;
            f_high = f_crit;
        } else {
            p_low = p_crit;
            f_low = f_crit;
        }
    }

    let g_crit = mass_flux_ps_eqm_throat(p_crit, s0);

    if debug {
        let region = ps_flash_region(p_crit, s0);
        let quality = x_ps_flash(p_crit, s0);
        dbg!(&(p_crit, region, quality, g_crit));
    }

    (p_crit, g_crit)
}
#[inline]
pub fn isentropic_pressure_scan_of_mass_flux(
    s0: SpecificHeatCapacity,
    p0: Pressure) -> () {

    let p_min_steam_table = Pressure::new::<megapascal>(0.000_611_212_677 * 1.01);
    
    // number of scan steps
    let n_steps = 1000;
    let dp = (p0 - p_min_steam_table) / n_steps as f64;

    let mut max_g_hem = MassFlux::ZERO;
    let mut p_crit = p0;

    let mut p_test = p0;

    for _ in 0..n_steps {
        p_test -= dp;

        if p_test < p_min_steam_table {
            break;
        }

        // skip region 1 — no choked flow in subcooled liquid
        let region = ps_flash_region(p_test, s0);
        if region == FwdEqnRegion::Region1 {
            continue;
        }

        let g_hem = mass_flux_ps_eqm_throat(p_test, s0);

        let quality = x_ps_flash(p_test, s0);
        dbg!(&(p_test, region, quality, g_hem, max_g_hem, p_crit));

        if g_hem > max_g_hem {
            max_g_hem = g_hem;
            p_crit = p_test;
            dbg!(&("new maximum found!", p_crit, max_g_hem));
        }
    }

}

/// estimates critical pressure ratio given ideal gas assumptions
/// for ideal gases, critical ratio depends on k 
/// but k is generally temperature dependent 
///
/// The evaluation here is to use throat properties to get the critical 
/// pressure ratio
///
#[inline]
pub fn get_critical_pressure_ratio_ideal_gas_using_throat_ph(
    p: Pressure,
    h: AvailableEnergy) -> Ratio {

    // note again that these are evaluated at throat
    let cp = cp_ph_eqm(p, h);
    let cv = cv_ph_eqm(p, h);

    let k = cp/cv;

    let ratio_one = Ratio::new::<ratio>(1.0);

    let k_plus_one = k + ratio_one;

    let k_minus_one = k - ratio_one;

    let exponent: f64 = (k/k_minus_one).get::<ratio>();
    let coeff: f64 = (2.0/k_plus_one).get::<ratio>();

    let ratio_value = coeff.powf(exponent);




    Ratio::new::<ratio>(ratio_value)
}


/// Finds the pressure where Mach number = 1 during isentropic expansion
/// This only works for superheated vapour
///
/// this takes in ph and optionally a stagnation entropy (if one wants to save 
/// on calculation speed)
#[inline]
pub fn get_critical_pressure_pure_vapour_ph_stagnation_properties(
    p0: Pressure,
    h0: AvailableEnergy,
    s0_opt: Option<SpecificHeatCapacity>) -> Pressure {

    let s0 = match s0_opt {
        Some(s0) => s0,
        None => s_ph_eqm(p0, h0),
    };

    // Initial guess: use ideal gas approximation as starting point
    //
    // Now, of course, the function below should use throat critical 
    // pressure, and I'm using stagnation properties 
    //
    // It is an approximation which works as a good initial guess 
    //
    // for moderate steam temperatures, 100-200C, cp/cv doesn't change much 
    // so that's fine. 
    //
    // For higher temperatures, this may not work as well 
    let p_guess = get_critical_pressure_ratio_ideal_gas_using_throat_ph(p0, h0)
        * p0;
    // ~(2/(k+1))^(k/(k-1)) for k≈1.3


    // Newton-Raphson or bisection to find where:
    // v = w (velocity equals speed of sound)
    //
    // From energy equation: h0 = h + v²/2
    // At critical point: v = w, so: h0 = h + w²/2

    let tolerance = Pressure::new::<pascal>(1.0); // 1 Pa tolerance
    let max_iterations = 50;

    // Bisection method bounds
    // Set bounds around the ideal gas guess (±30% to be safe)
    // This reduces iterations compared to starting at 0.1*p0 to 1.0*p0
    let mut p_low = p_guess * 0.7;   // 30% below guess
    let mut p_high = p_guess * 1.3;  // 30% above guess

    // Clamp bounds to reasonable range
    if p_low < p0 * 0.1 {
        p_low = p0 * 0.1;
    }
    if p_high > p0 * 0.99 {
        p_high = p0 * 0.99;
    }


    for _ in 0..max_iterations {
        let p_mid = (p_low + p_high) / 2.0;

        // Get properties at this pressure (isentropic)
        let h_mid = h_ps_eqm(p_mid, s0);
        let w_mid = w_ph_wood_wallis(p_mid, h_mid);

        // Calculate velocity from energy equation
        // h0 = h + v²/2  =>  v = sqrt(2*(h0 - h))
        let delta_h = h0 - h_mid;

        if delta_h < AvailableEnergy::ZERO {
            // Pressure too low, expansion exceeded stagnation enthalpy
            p_low = p_mid;
            continue;
        }

        let v_squared = 2.0 * delta_h;
        let v = v_squared.sqrt();

        // Check if Mach = 1 (v = w)
        let mach = v / w_mid;
        let mach_value = mach.get::<ratio>();

        if (mach_value - 1.0).abs() < 0.0001 {
            return p_mid;
        }

        // Adjust bounds
        if mach_value < 1.0 {
            p_high = p_mid; // Need lower pressure (more expansion)
        } else {
            p_low = p_mid;  // Need higher pressure (less expansion)
        }

        // Check convergence
        if (p_high - p_low) < tolerance {
            return (p_low + p_high) / 2.0;
        }
    }

    // Return midpoint if not converged
    (p_low + p_high) / 2.0
}
