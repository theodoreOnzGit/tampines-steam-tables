use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, SystemTime};

use teh_o_prke::decay_heat::DecayHeat;
use teh_o_prke::feedback_mechanisms::fission_product_poisons::Xenon135Poisoning;
use teh_o_prke::zero_power_prke::six_group_precursor_prke::six_group_constants::FissioningNuclideType;
use teh_o_prke::{feedback_mechanisms::SixFactorFormulaFeedback, zero_power_prke::six_group_precursor_prke::SixGroupPRKE};
use uom::si::area::square_meter;
use uom::si::energy::{kilojoule, megaelectronvolt};
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::linear_number_density::per_meter;
use uom::si::mass::kilogram;
use uom::si::power::megawatt;
use uom::si::time::{microsecond, second};
use uom::si::velocity::meter_per_second;
use uom::si::volume::cubic_meter;
use uom::si::volumetric_number_rate::per_cubic_meter_second;
use uom::si::{f64::*, ratio::ratio};
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::ConstZero;

use crate::{FHRSimulatorApp, FHRState};

impl FHRSimulatorApp {

    /// associated function for PRKE calculation 
    /// for continuous loop
    pub fn calculate_prke_loop(
        fhr_state: Arc<Mutex<FHRState>>){

        // construct a prke six group object
        // probably want to use a u235 group or u233 group
        // default
        let mut prke_six_group :SixGroupPRKE = SixGroupPRKE::default();

        let prke_timestep = Time::new::<microsecond>(25.0);
        let reactor_volume = Volume::new::<cubic_meter>(0.5);
        let macroscopic_fission_xs = LinearNumberDensity::new::<per_meter>(1.0);
        let mut pebble_bed_th_struct = 
            PebbleBedThermalHydraulics::new();
        let fhr_state_clone = fhr_state.clone();

        // then decay heat struct 
        let mut fhr_decay_heat_struct = DecayHeat::default();

        // then xenon poisoning struct 
        let mut fhr_xe135_poisoning = Xenon135Poisoning::default();


        // now, time controls 
        let loop_time = SystemTime::now();
        let mut current_simulation_time = Time::ZERO;

        // calculation loop (indefinite)
        //
        //// calculation time and time to sleep
        // to be done once every timestep
        loop {

            // so now, let's do the necessary things
            // first, timestep and loop time 
            //
            // second, read and update the local_ciet_state

            let loop_time_start = loop_time.elapsed().unwrap();
            // now this is arbitrary, user can set
            let mut keff_six_factor = SixFactorFormulaFeedback::default();
            // start from keff = 1
            // all terms = 1
            // start with leakage
            keff_six_factor.p_tnl = Ratio::new::<ratio>(0.9);
            keff_six_factor.p_fnl = Ratio::new::<ratio>(0.7);
            // then fuel reproduction
            keff_six_factor.eta = Ratio::new::<ratio>(2.2);
            // resonance esc probability 
            keff_six_factor.p = Ratio::new::<ratio>(0.8);
            // thermal utilisation 
            keff_six_factor.f = Ratio::new::<ratio>(0.9);

            
            // fast fission 
            keff_six_factor.epsilon = Ratio::new::<ratio>(1.03);
            // keff total is about 1.0278
            // excess reactivity is about 0.0278
            //
            // basically control rod should be about this much 
            // and fuel temp feedback about this much also
            // some of these are arbitrary
            Self::calculate_prke_for_one_timestep(
                &mut fhr_state_clone.lock().unwrap(),
                &mut keff_six_factor,
                &mut prke_six_group,
                prke_timestep,
                reactor_volume,
                macroscopic_fission_xs,
                &mut fhr_decay_heat_struct,
                &mut pebble_bed_th_struct,
                &mut fhr_xe135_poisoning,
            );



            current_simulation_time += prke_timestep;

            let prke_simulation_time_seconds = current_simulation_time.get::<second>();

            let prke_elapsed_time_seconds = 
                (loop_time.elapsed().unwrap().as_secs_f64() * 100.0).round()/100.0;

            let overall_simulation_in_realtime_or_faster: bool = 
                prke_simulation_time_seconds >= prke_elapsed_time_seconds;

            // now update the fhr state 
            let prke_timestep_microseconds = prke_timestep.get::<microsecond>();

            fhr_state_clone.lock().unwrap().prke_timestep_microseconds 
                = prke_timestep_microseconds;





            fhr_state_clone.lock().unwrap().prke_simulation_time_seconds 
                = prke_simulation_time_seconds;

            fhr_state_clone.lock().unwrap().prke_elapsed_time_seconds 
                = prke_elapsed_time_seconds;


            // calculation time and time to sleep
            let loop_time_end = loop_time.elapsed().unwrap();
            let time_taken_for_calculation_loop_microseconds: f64 = 
                (loop_time_end - loop_time_start)
                .as_micros() as f64;
            fhr_state_clone.lock().unwrap().prke_calc_time_microseconds 
                = time_taken_for_calculation_loop_microseconds;

            let time_to_sleep_microseconds: u64 = 
                (prke_timestep.get::<microsecond>() - 
                 time_taken_for_calculation_loop_microseconds)
                .round().abs() as u64;

            let time_to_sleep: Duration = 
                Duration::from_micros(time_to_sleep_microseconds - 1);


            // last condition for sleeping
            let real_time_in_current_timestep: bool = 
                time_to_sleep_microseconds > 1;

            //
            let fast_forward_botton_on = false;

            if overall_simulation_in_realtime_or_faster && 
                real_time_in_current_timestep && 
                    !fast_forward_botton_on 
            {
                thread::sleep(time_to_sleep);
            } else if overall_simulation_in_realtime_or_faster 
                && real_time_in_current_timestep 
                    && fast_forward_botton_on 
            {
                // sleep 5 microseconds if fast fwd
                let short_time_to_sleep: Duration = Duration::from_micros(5);
                thread::sleep(short_time_to_sleep);
            } else {
                // don't sleep

            }
            //let time_to_sleep = Duration::from_millis(40);

            //dbg!(&(
            //        time_taken_for_calculation_loop_microseconds,
            //        prke_timestep.get::<microsecond>(),
            //)
            //);
            //thread::sleep(time_to_sleep);
        }

    }
    /// associated function for PRKE calculation 
    /// for single timestep
    /// for prke anyway
    /// note that the prke timestep will be different from 
    /// the main thermal hydraulics timestep
    pub fn calculate_prke_for_one_timestep(
        fhr_state_ref: &mut FHRState,
        keff_six_factor: &mut SixFactorFormulaFeedback,
        prke_six_group: &mut SixGroupPRKE,
        prke_timestep: Time,
        reactor_volume: Volume,
        macroscopic_fission_xs: LinearNumberDensity,
        fhr_decay_heat: &mut DecayHeat,
        pebble_bed_th_struct: &mut PebbleBedThermalHydraulics,
        fhr_xe135_poisoning: &mut Xenon135Poisoning,
        ){

        // within each timestep, I need to obtain feedback
        // so basically for now, fuel temperature and control rod insertion 
        // based on that, calculate PRKE


        let fuel_temp: ThermodynamicTemperature = 
            ThermodynamicTemperature::new::<degree_celsius>(
                fhr_state_ref.pebble_core_temp_degc
            );
        let left_cr_insertion_frac = 
            fhr_state_ref.left_cr_insertion_frac;
        let right_cr_insertion_frac = 
            fhr_state_ref.right_cr_insertion_frac;

        // now based on this, calculate feedback
        //
        let left_cr_insertion_ratio = 
            Ratio::new::<ratio>(left_cr_insertion_frac as f64);
        let right_cr_insertion_ratio = 
            Ratio::new::<ratio>(right_cr_insertion_frac as f64);

        keff_six_factor.fuel_temp_feedback(
            fuel_temp, 
            FHRSimulatorApp::fuel_temp_resonance_esc_feedback_linear
        );
        keff_six_factor.control_rod_feedback(
            left_cr_insertion_ratio, 
            FHRSimulatorApp::fuel_utilisation_factor_chg_for_control_rod_polynomial
        );
        keff_six_factor.control_rod_feedback(
            right_cr_insertion_ratio, 
            FHRSimulatorApp::fuel_utilisation_factor_chg_for_control_rod_polynomial
        );

        // next xenon poisoning feedback
        //
        //
        // adjust for xenon poisoning
        let xe135_mass_conc = fhr_xe135_poisoning.get_current_xe135_conc();
        let thermal_utilisation_feedback_fractional_chg_from_xenon: f64 = 
            Xenon135Poisoning::simplified_poison_concentration_feedback(
                xe135_mass_conc
            ).get::<ratio>();

        keff_six_factor.f *= 
            thermal_utilisation_feedback_fractional_chg_from_xenon;


        // after feedback we should get the reactivity 
        let reactivity: Ratio = keff_six_factor.calc_rho();
        let neutron_generation_time = Time::new::<second>(2.31e-4);
        let mean_neutron_time = neutron_generation_time/keff_six_factor.calc_keff();
        let background_source_rate = 
            VolumetricNumberRate::new::<per_cubic_meter_second>(5.0);

        // use this to toggle between implicit and explicit solver
        let implicit_solver = false;


        if implicit_solver == true {
            let _neutron_pop_and_six_group_precursor_vec = 
                prke_six_group.solve_next_timestep_precursor_concentration_and_neutron_pop_vector_implicit(
                    prke_timestep, 
                    reactivity, 
                    mean_neutron_time, 
                    background_source_rate
                );
        } else {

            let _neutron_pop_and_six_group_precursor_vec = 
                prke_six_group.solve_next_timestep_precursor_concentration_and_neutron_pop_vector_explicit(
                    prke_timestep, 
                    reactivity, 
                    mean_neutron_time, 
                    background_source_rate
                );
        }

        // then get the current neutron population 
        let current_neutron_pop_density: VolumetricNumberDensity = 
            prke_six_group.get_current_neutron_population_density();
        // then total reactor volume to get neutron number 
        //let current_neutron_pop: Ratio =  
        //    reactor_volume * current_neutron_pop_density;
        

        // we can get a power production 
        // and some of it should go to decay heat
        // so we need Sigma_f * phi 
        let neutron_speed: Velocity = Velocity::new::<meter_per_second>(2200.0);
        let current_neutron_flux: ArealNumberRate = 
            (current_neutron_pop_density * neutron_speed).into();

        // then fission rate
        // should be a number Rate
        // per unit vol
        let fission_rate_density: VolumetricNumberRate = 
            (current_neutron_flux * macroscopic_fission_xs).into();

        // should be a frequency
        let fission_rate: Frequency  = 
            fission_rate_density * reactor_volume;

        // note: it's convenient here to calc xe135 feedback 

        let fissioning_nuclide = FissioningNuclideType::U235;
        let _xe135_conc_next_timestep = 
            fhr_xe135_poisoning.calc_xe_135_and_return_num_density(
                prke_timestep, 
                fission_rate_density, 
                fissioning_nuclide, 
                current_neutron_pop_density);

        let power_per_fission = 
            Energy::new::<megaelectronvolt>(200.0);

        // immediate power from fission
        let fission_power_instantaneous: Power = 
            power_per_fission * fission_rate;

        // add to decay heat precursors 
        fhr_decay_heat.add_decay_heat_precursor1(
            fission_power_instantaneous * 0.04, prke_timestep
        );
        fhr_decay_heat.add_decay_heat_precursor2(
            fission_power_instantaneous * 0.04, prke_timestep
        );
        fhr_decay_heat.add_decay_heat_precursor3(
            fission_power_instantaneous * 0.02, prke_timestep
        );


        // adjust fission power for decay heat 
        // fission power less decay heat = 1.0 - 0.04 - 0.04 - 0.02 = 0.9
        let mut fission_power_corrected_for_decay_heat = fission_power_instantaneous * 0.9;
        let mut reactor_current_decay_heat: Power = fhr_decay_heat.calc_decay_heat_power_1(prke_timestep).abs();
        reactor_current_decay_heat += fhr_decay_heat.calc_decay_heat_power_2(prke_timestep).abs();
        reactor_current_decay_heat += fhr_decay_heat.calc_decay_heat_power_3(prke_timestep).abs();
        fission_power_corrected_for_decay_heat += reactor_current_decay_heat;

        // with the correct fission power now, we can 
        // calc temperature
        //
        //
        // These are arbitrary values, will adjust later

        let pebble_bed_mass = Mass::new::<kilogram>(8000.0);
        let pebble_bed_heat_transfer_area = Area::new::<square_meter>(300.0);
        let pebble_bed_overall_htc = HeatTransfer::new::<watt_per_square_meter_kelvin>(400.0);
        let pebble_bed_coolant_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_ref.pebble_bed_coolant_temp_degc
        );

        let heat_removal_from_pebble_bed = 
            pebble_bed_th_struct.calc_th_and_return_heat_removal_from_pebble_bed(
                prke_timestep, 
                fission_power_corrected_for_decay_heat, 
                pebble_bed_mass, 
                pebble_bed_heat_transfer_area, 
                pebble_bed_overall_htc, 
                pebble_bed_coolant_temp);

        let pebble_bed_fuel_temp = pebble_bed_th_struct.get_temperature_from_enthalpy_uo2_heuristic();

        // update the fhr state 
        fhr_state_ref.pebble_core_temp_degc = 
            pebble_bed_fuel_temp.get::<degree_celsius>();

        // the heat removal from pebble bed should be stored in fhr 
        // state, so as to sync correctly with the coolant
        // the prke timestep will also be added so as to obtain an 
        // average heat removal rate to be transferred to the coolant

        fhr_state_ref.prke_loop_accumulated_heat_removal_kilojoules
            += (heat_removal_from_pebble_bed*prke_timestep).get::<kilojoule>();
        fhr_state_ref.prke_loop_accumulated_timestep_seconds
            += prke_timestep.get::<second>();

        // reactor power 
        //

        let keff = keff_six_factor.calc_keff();
        fhr_state_ref.keff = keff.get::<ratio>();
        fhr_state_ref.reactor_power_megawatts = 
            fission_power_corrected_for_decay_heat.get::<megawatt>();
        fhr_state_ref.reactor_decay_heat_megawatts = 
            reactor_current_decay_heat.get::<megawatt>();

        // reactivity in dollars 
        let beta_delayed_frac_total = prke_six_group.get_total_delayed_fraction();
        let reactivity_dollars: f64 
            = (reactivity/beta_delayed_frac_total).get::<ratio>();

        fhr_state_ref.reactivity_dollars = reactivity_dollars;

        let xenon135_feedback_dollars_approx = 
            (thermal_utilisation_feedback_fractional_chg_from_xenon-1.0)/
            thermal_utilisation_feedback_fractional_chg_from_xenon/
            beta_delayed_frac_total;

        fhr_state_ref.xenon135_feedback_dollars = 
            xenon135_feedback_dollars_approx.get::<ratio>();


        let debug_settings = false;
        if debug_settings {
            // that settles thermal hydraulics
            dbg!(&(
                    pebble_bed_fuel_temp,
                    keff,
                    reactivity,
                    fission_power_instantaneous,
                    fission_power_corrected_for_decay_heat,
                    heat_removal_from_pebble_bed
            ));
        }


    }

    /// fuel temperature feedback and how it deals with resonance escape 
    /// probability 
    ///
    /// basically, with increased fuel temp, lower resonance 
    /// esc probability
    ///
    /// this is a heuristic, can be replaced by more complex 
    /// functions later
    ///
    /// the ratio u return is factor of increase or decrease in resonance 
    /// escape probability
    pub fn fuel_temp_resonance_esc_feedback_linear(fuel_temp: ThermodynamicTemperature) -> Ratio {

        let fuel_temp_degc = fuel_temp.get::<degree_celsius>();
        // from an arbitrary map of:
        // fuel_temp_degc, resonance esc change factor
        // 200,1.2,
        // 300,1.1,
        // 400,1.05,
        // 500,1,
        // 600,0.99,
        // 700,0.98,
        // 800,0.97,
        // 900,0.96,
        // 1000,0.95,
        // 1100,0.94,
        // 1200,0.9,
        // 1300,0.85,
        // 1400,0.83,
        // 1500,0.8,
        //
        // I get a trendline in LibreOffice

        return Ratio::new::<ratio>(-2.3077e-4 * fuel_temp_degc + 1.2462e0);


    }

    pub fn fuel_utilisation_factor_chg_for_control_rod_polynomial(
        cr_insertion_factor: Ratio) -> Ratio {

        //this is an arbitrary map, but just useful for 
        //simulation
        //control rod insertion frac,fuel utilisation factor change
        // 0,1.03
        // 0.1,1.0275
        // 0.2,1.02
        // 0.3,1.005
        // 0.4,0.975
        // 0.5,0.875
        // 0.6,0.82
        // 0.7,0.775
        // 0.8,0.77
        // 0.9,0.76
        // 1,0.75

        let cr_insertion_factor_f64: f64 = 
            cr_insertion_factor.get::<ratio>();
        let term1 = -8.413e0 * cr_insertion_factor_f64.powi(5);
        let term2 = 2.064e1 * cr_insertion_factor_f64.powi(4);
        let term3 = -1.639e1 * cr_insertion_factor_f64.powi(3);
        let term4 = 4.238e0 * cr_insertion_factor_f64.powi(2);
        let term5 = -3.626e-1 * cr_insertion_factor_f64.powi(1);
        let term6 = 1.031e0 * cr_insertion_factor_f64.powi(0);


        return Ratio::new::<ratio>(term1 + term2 + term3 + term4 + term5 + term6);
    }
    
    // this is a dummy function so that no 
    // matter what you input, no change is given
    pub fn constant_ratio_no_change_function(_: Ratio) -> Ratio {
        return Ratio::new::<ratio>(1.0);
    }
    
}



/// this contains code for TIGHTLY coupled thermal hydraulics,
/// that is between the PRKE
pub mod pebble_bed_thermal_hydraulics;
pub use pebble_bed_thermal_hydraulics::*;


