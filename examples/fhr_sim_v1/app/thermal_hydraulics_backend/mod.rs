use std::sync::{Arc, Mutex};

use std::thread;
use std::time::{Duration, SystemTime};

use fhr_thermal_hydraulics_state::FHRThermalHydraulicsState;
use ndarray::{Array, Array1};
use tuas_boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
use tuas_boussinesq_solver::pre_built_components::shell_and_tube_heat_exchanger::SimpleShellAndTubeHeatExchanger;
use tuas_boussinesq_solver::prelude::beta_testing::{FluidArray, HeatTransferEntity, HeatTransferInteractionType};
use uom::si::energy::kilojoule;
use uom::si::mass_rate::kilogram_per_second;
use uom::si::power::megawatt;
use uom::si::pressure::{kilopascal, megapascal};
use uom::si::thermal_conductance::watt_per_kelvin;
//use teh_o_prke::decay_heat::DecayHeat;
//use teh_o_prke::feedback_mechanisms::fission_product_poisons::Xenon135Poisoning;
//use teh_o_prke::zero_power_prke::six_group::FissioningNuclideType;
//use teh_o_prke::{feedback_mechanisms::SixFactorFormulaFeedback, zero_power_prke::six_group::SixGroupPRKE};
//use uom::si::area::square_meter;
//use uom::si::energy::{kilojoule, megaelectronvolt};
//use uom::si::heat_transfer::watt_per_square_meter_kelvin;
//use uom::si::linear_number_density::per_meter;
//use uom::si::mass::kilogram;
//use uom::si::power::megawatt;
use uom::si::time::{microsecond, second};
//use uom::si::velocity::meter_per_second;
//use uom::si::volume::cubic_meter;
//use uom::si::volumetric_number_rate::per_cubic_meter_second;
use uom::si::f64::*;
//use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::ConstZero;


use components::*;
use pri_loop_fluid_mechanics_calc_fns::four_branch_pri_and_intermediate_loop_fluid_mechanics_only;
use tuas_boussinesq_solver::pre_built_components::insulated_pipes_and_fluid_components::InsulatedFluidComponent;
use tuas_boussinesq_solver::pre_built_components::non_insulated_fluid_components::NonInsulatedFluidComponent;
use crate::{FHRSimulatorApp, FHRState};



impl FHRSimulatorApp {

    /// for the gFHR primary loop, and intermediate loop 
    /// there are four branches that need to be solved for flowrate 
    ///
    /// this code handles the solution procedure
    /// using the tuas_boussinesq_solver library code
    ///
    /// and handles fluid mechanics and heat transfer for one time step
    pub(crate) fn four_branch_pri_and_intermediate_loop_single_time_step(
        pri_loop_pump_pressure: Pressure,
        intrmd_loop_pump_pressure: Pressure,
        reactor_power: Power,
        timestep: Time,
        // diagnostics 
        simulation_time: Time,
        // reactor branch
        reactor_pipe_1: &mut InsulatedFluidComponent,
        // downcomer branch 1
        downcomer_pipe_2: &mut InsulatedFluidComponent,
        // downcomer branch 2
        downcomer_pipe_3: &mut InsulatedFluidComponent,
        // mixing nodes for pri loop
        bottom_mixing_node_pri_loop: &mut HeatTransferEntity,
        top_mixing_node_pri_loop: &mut HeatTransferEntity,
        // Intermediate heat exchanger branch in pri loop
        fhr_pipe_11: &mut InsulatedFluidComponent,
        fhr_pipe_10: &mut InsulatedFluidComponent,
        fhr_pri_loop_pump_9: &mut NonInsulatedFluidComponent,
        fhr_pipe_8: &mut InsulatedFluidComponent,
        fhr_pipe_7: &mut InsulatedFluidComponent,
        ihx_sthe_6: &mut SimpleShellAndTubeHeatExchanger,
        fhr_pipe_5: &mut InsulatedFluidComponent,
        fhr_pipe_4: &mut InsulatedFluidComponent,
        // intermediate loop ihx side
        fhr_pipe_17: &mut InsulatedFluidComponent,
        fhr_pipe_12: &mut InsulatedFluidComponent,
        // intermediate loop steam generator side
        fhr_intrmd_loop_pump_16: &mut NonInsulatedFluidComponent,
        fhr_pipe_15: &mut InsulatedFluidComponent,
        fhr_steam_generator_shell_side_14: &mut NonInsulatedFluidComponent,
        fhr_pipe_13: &mut InsulatedFluidComponent,
        // mixing nodes for intermediate loop
        bottom_mixing_node_intrmd_loop: &mut HeatTransferEntity,
        top_mixing_node_intrmd_loop: &mut HeatTransferEntity,
        // steam generator settings 
        steam_generator_tube_side_temperature: ThermodynamicTemperature,
        steam_generator_overall_ua: ThermalConductance,

        ) -> FHRThermalHydraulicsState {

            // fluid mechnaics portion for both loops


            let (reactor_branch_flow, downcomer_branch_1_flow,
                downcomer_branch_2_flow, pri_loop_intermediate_heat_exchanger_branch_flow,
                intrmd_loop_ihx_br_flow, intrmd_loop_steam_gen_br_flow)
                = four_branch_pri_and_intermediate_loop_fluid_mechanics_only(
                    pri_loop_pump_pressure, 
                    intrmd_loop_pump_pressure, 
                    reactor_pipe_1, 
                    downcomer_pipe_2, 
                    downcomer_pipe_3, 
                    fhr_pipe_11, 
                    fhr_pipe_10, 
                    fhr_pri_loop_pump_9, 
                    fhr_pipe_8, 
                    fhr_pipe_7, 
                    ihx_sthe_6, 
                    fhr_pipe_5, 
                    fhr_pipe_4, 
                    fhr_pipe_17, 
                    fhr_pipe_12, 
                    fhr_intrmd_loop_pump_16, 
                    fhr_pipe_15, 
                    fhr_steam_generator_shell_side_14, 
                    fhr_pipe_13);

            // thermal hydraulics part
            //
            // first, we are going to make heat transfer interactions
            // need rough temperature for density calcs, not super 
            // important as we assume boussineq approximation 
            // ie density differences only important for buoyancy calcs
            let average_temperature_for_density_calcs_pri_loop = 
                ThermodynamicTemperature::new::<degree_celsius>(600.0);


            let average_flibe_density = 
                LiquidMaterial::FLiBe.try_get_density(
                    average_temperature_for_density_calcs_pri_loop).unwrap();

            let downcomer_branch_1_advection_heat_transfer_interaction = 
                HeatTransferInteractionType::
                new_advection_interaction(downcomer_branch_1_flow, 
                    average_flibe_density, 
                    average_flibe_density);

            let downcomer_branch_2_advection_heat_transfer_interaction = 
                HeatTransferInteractionType::
                new_advection_interaction(downcomer_branch_2_flow, 
                    average_flibe_density, 
                    average_flibe_density);

            let reactor_branch_advection_heat_transfer_interaction = 
                HeatTransferInteractionType::
                new_advection_interaction(reactor_branch_flow, 
                    average_flibe_density, 
                    average_flibe_density);

            let ihx_advection_heat_transfer_interaction = 
                HeatTransferInteractionType::
                new_advection_interaction(pri_loop_intermediate_heat_exchanger_branch_flow, 
                    average_flibe_density, 
                    average_flibe_density);
            // for intermediate loop, we use lower temp, 
            // about 450 C
            //
            // as it is a HITEC salt (nitrate salt)
            let average_temperature_for_density_calcs_intrmd_loop = 
                ThermodynamicTemperature::new::<degree_celsius>(450.0);

            let average_hitec_density = 
                LiquidMaterial::HITEC.try_get_density(
                    average_temperature_for_density_calcs_intrmd_loop).unwrap();

            let intrmd_loop_ihx_br_heat_transfer_interaction = 
                HeatTransferInteractionType::
                new_advection_interaction(intrmd_loop_ihx_br_flow, 
                    average_hitec_density, 
                    average_hitec_density);
            let intrmd_loop_steam_gen_br_heat_transfer_interaction = 
                HeatTransferInteractionType::
                new_advection_interaction(intrmd_loop_steam_gen_br_flow, 
                    average_hitec_density, 
                    average_hitec_density);

            // note that reactor branch flow, 
            // downcomer_branch_1_flow, 
            // downcomer_branch_1_flow and 
            // intermediate_heat_exchanger_branch_flow in the pri loop 
            // all go from bottom mixing node to top mixing node
            //
            // with this in mind, we now link up the components 

            // downcomer 1 branch
            {
                bottom_mixing_node_pri_loop.link_to_front(
                    &mut downcomer_pipe_2.pipe_fluid_array, 
                    downcomer_branch_1_advection_heat_transfer_interaction)
                    .unwrap();

                downcomer_pipe_2.pipe_fluid_array.link_to_front(
                    top_mixing_node_pri_loop, 
                    downcomer_branch_1_advection_heat_transfer_interaction)
                    .unwrap();

                }
            // downcomer 2 branch
            {
                bottom_mixing_node_pri_loop.link_to_front(
                    &mut downcomer_pipe_3.pipe_fluid_array, 
                    downcomer_branch_2_advection_heat_transfer_interaction)
                    .unwrap();

                downcomer_pipe_3.pipe_fluid_array.link_to_front(
                    top_mixing_node_pri_loop, 
                    downcomer_branch_2_advection_heat_transfer_interaction)
                    .unwrap();
                }
            // pri loop 
            // ihx branch 
            {

                bottom_mixing_node_pri_loop.link_to_front(
                    &mut fhr_pipe_11.pipe_fluid_array, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_11.pipe_fluid_array.link_to_front(
                    &mut fhr_pipe_10.pipe_fluid_array, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_10.pipe_fluid_array.link_to_front(
                    &mut fhr_pri_loop_pump_9.pipe_fluid_array, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pri_loop_pump_9.pipe_fluid_array.link_to_front(
                    &mut fhr_pipe_8.pipe_fluid_array, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_8.pipe_fluid_array.link_to_front(
                    &mut fhr_pipe_7.pipe_fluid_array, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_7.pipe_fluid_array.link_to_front(
                    &mut ihx_sthe_6.shell_side_fluid_array, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                ihx_sthe_6.shell_side_fluid_array.link_to_front(
                    &mut fhr_pipe_5.pipe_fluid_array,
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_5.pipe_fluid_array.link_to_front(
                    &mut fhr_pipe_4.pipe_fluid_array,
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_4.pipe_fluid_array.link_to_front(
                    top_mixing_node_pri_loop, 
                    ihx_advection_heat_transfer_interaction)
                    .unwrap();
                }

            // intermediate loop ihx branch 
            {

                bottom_mixing_node_intrmd_loop.link_to_front(
                    &mut fhr_pipe_17.pipe_fluid_array, 
                    intrmd_loop_ihx_br_heat_transfer_interaction)
                    .unwrap();

                ihx_sthe_6.tube_side_fluid_array_for_single_tube.link_to_back(
                    &mut fhr_pipe_17.pipe_fluid_array, 
                    intrmd_loop_ihx_br_heat_transfer_interaction)
                    .unwrap();

                ihx_sthe_6.tube_side_fluid_array_for_single_tube.link_to_front(
                    &mut fhr_pipe_12.pipe_fluid_array, 
                    intrmd_loop_ihx_br_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_12.pipe_fluid_array.link_to_front(
                    top_mixing_node_intrmd_loop, 
                    intrmd_loop_ihx_br_heat_transfer_interaction)
                    .unwrap();

                }

            // intermediate loop steam generator branch
            {

                bottom_mixing_node_intrmd_loop.link_to_front(
                    &mut fhr_intrmd_loop_pump_16.pipe_fluid_array, 
                    intrmd_loop_steam_gen_br_heat_transfer_interaction)
                    .unwrap();

                fhr_intrmd_loop_pump_16.pipe_fluid_array.link_to_front(
                    &mut fhr_pipe_15.pipe_fluid_array, 
                    intrmd_loop_steam_gen_br_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_15.pipe_fluid_array.link_to_front(
                    &mut fhr_steam_generator_shell_side_14.pipe_fluid_array, 
                    intrmd_loop_steam_gen_br_heat_transfer_interaction)
                    .unwrap();

                fhr_steam_generator_shell_side_14.pipe_fluid_array.link_to_front(
                    &mut fhr_pipe_13.pipe_fluid_array, 
                    intrmd_loop_steam_gen_br_heat_transfer_interaction)
                    .unwrap();

                fhr_pipe_13.pipe_fluid_array.link_to_front(
                    top_mixing_node_intrmd_loop, 
                    intrmd_loop_steam_gen_br_heat_transfer_interaction)
                    .unwrap();
                }
            {

                // for steam generator, I want to manually remove heat from it 
                // uniformly 

                let number_of_temperature_nodes_for_sg = 2;
                let mut q_frac_arr: Array1<f64> = Array::default(number_of_temperature_nodes_for_sg);
                // we want the middle node to contain all the power
                q_frac_arr[0] = 0.5;
                q_frac_arr[1] = 0.5;
                let mut sg_fluid_array_clone: FluidArray = 
                    fhr_steam_generator_shell_side_14
                    .pipe_fluid_array
                    .clone()
                    .try_into()
                    .unwrap();
                let steam_gen_heat_change: Power;

                let temperature_diff = 
                    TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(
                        sg_fluid_array_clone.try_get_bulk_temperature()
                        .unwrap()
                        .get::<degree_celsius>() 
                        - steam_generator_tube_side_temperature
                        .get::<degree_celsius>()
                    );

                // Q_added_to_destination = -UA*(T_destination - T_steam)
                steam_gen_heat_change = -temperature_diff*steam_generator_overall_ua;

                sg_fluid_array_clone
                    .lateral_link_new_power_vector(
                        steam_gen_heat_change, 
                        q_frac_arr)
                    .unwrap();

                fhr_steam_generator_shell_side_14.pipe_fluid_array
                    = sg_fluid_array_clone.into();
            }
            // now for the reactor branch, we must have some kind of 
            // power input here 
            {

                // i'll use the lateral link new power vector code 
                //
                // this sets the reactor power in the middle part of the 
                // pipe
                let number_of_temperature_nodes_for_reactor = 5;
                let mut q_frac_arr: Array1<f64> = Array::default(number_of_temperature_nodes_for_reactor);
                // we want the middle node to contain all the power
                q_frac_arr[0] = 0.0;
                q_frac_arr[1] = 0.0;
                q_frac_arr[2] = 1.0;
                q_frac_arr[3] = 0.0;
                q_frac_arr[4] = 0.0;

                // now i need to get the fluid array out first 

                let mut reactor_fluid_array_clone: FluidArray = 
                    reactor_pipe_1
                    .pipe_fluid_array
                    .clone()
                    .try_into()
                    .unwrap();

                reactor_fluid_array_clone
                    .lateral_link_new_power_vector(
                        reactor_power, 
                        q_frac_arr)
                    .unwrap();

                reactor_pipe_1.pipe_fluid_array = 
                    reactor_fluid_array_clone.into();

                // now, add the connections

                reactor_pipe_1.pipe_fluid_array.link_to_front(
                    top_mixing_node_pri_loop, 
                    reactor_branch_advection_heat_transfer_interaction)
                    .unwrap();
                reactor_pipe_1.pipe_fluid_array.link_to_back(
                    bottom_mixing_node_pri_loop, 
                    reactor_branch_advection_heat_transfer_interaction)
                    .unwrap();
                }

            // now we are ready to advance timesteps for all components 
            // and mixing nodes 

            let zero_power = Power::ZERO;
            // for pri loop 
            // I'm not going to add another round of power 
            // because I already added it to the top
            // so i'll just add zero power
            //
            // this is reactor and downcomer branches
            {
                reactor_pipe_1
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        reactor_branch_flow, 
                        zero_power)
                    .unwrap();

                downcomer_pipe_2
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        downcomer_branch_1_flow, 
                        zero_power)
                    .unwrap();

                downcomer_pipe_3
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        downcomer_branch_2_flow, 
                        zero_power)
                    .unwrap();
                }

            // this is the pri loop ihx branch
            // except for the ihx itself
            {

                fhr_pipe_11
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_10
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                fhr_pri_loop_pump_9
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_8
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_7
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_5
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_4
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        pri_loop_intermediate_heat_exchanger_branch_flow, 
                        zero_power)
                    .unwrap();
                }

            // ihx 
            {

                let prandtl_wall_correction_setting = true; 
                let tube_side_total_mass_flowrate = intrmd_loop_ihx_br_flow;
                let shell_side_total_mass_flowrate = pri_loop_intermediate_heat_exchanger_branch_flow;

                ihx_sthe_6.lateral_and_miscellaneous_connections(
                    prandtl_wall_correction_setting, 
                    tube_side_total_mass_flowrate, 
                    shell_side_total_mass_flowrate).unwrap();

            }
            // hitec intrmd loop 
            //
            // except for ihx itself
            {
                // ihx branch
                fhr_pipe_17
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        intrmd_loop_ihx_br_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_12
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        intrmd_loop_ihx_br_flow, 
                        zero_power)
                    .unwrap();

                // steam gen branch
                fhr_intrmd_loop_pump_16
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        intrmd_loop_steam_gen_br_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_15
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        intrmd_loop_steam_gen_br_flow, 
                        zero_power)
                    .unwrap();
                fhr_steam_generator_shell_side_14
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        intrmd_loop_steam_gen_br_flow, 
                        zero_power)
                    .unwrap();
                fhr_pipe_13
                    .lateral_and_miscellaneous_connections_no_wall_correction(
                        intrmd_loop_steam_gen_br_flow, 
                        zero_power)
                    .unwrap();
                }

            // timestep advance for all heat transfer entities
            {
                // pri loop (with ihx)
                reactor_pipe_1
                    .advance_timestep(timestep)
                    .unwrap();
                downcomer_pipe_2
                    .advance_timestep(timestep)
                    .unwrap();
                downcomer_pipe_3
                    .advance_timestep(timestep)
                    .unwrap();


                fhr_pipe_4
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_5
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_7
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_8
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pri_loop_pump_9
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_10
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_11
                    .advance_timestep(timestep)
                    .unwrap();

                // intermediate branch (less ihx)
                fhr_pipe_12
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_17
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_13
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_steam_generator_shell_side_14
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_pipe_15
                    .advance_timestep(timestep)
                    .unwrap();
                fhr_intrmd_loop_pump_16
                    .advance_timestep(timestep)
                    .unwrap();

                // all mixing nodes
                top_mixing_node_pri_loop
                    .advance_timestep_mut_self(timestep)
                    .unwrap();
                bottom_mixing_node_pri_loop
                    .advance_timestep_mut_self(timestep)
                    .unwrap();
                top_mixing_node_intrmd_loop
                    .advance_timestep_mut_self(timestep)
                    .unwrap();
                bottom_mixing_node_intrmd_loop
                    .advance_timestep_mut_self(timestep)
                    .unwrap();

                ihx_sthe_6
                    .advance_timestep(timestep)
                    .unwrap();
                }

            // now I want reactor temperature profile 
            let reactor_temp_profile: Vec<ThermodynamicTemperature> = 
                reactor_pipe_1
                .pipe_fluid_array_temperature()
                .unwrap();
            let reactor_temp_profile_degc: Vec<f64> = 
                reactor_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // sthe temperature profile
            let ihx_shell_side_temp_profile: Vec<ThermodynamicTemperature> = 
                ihx_sthe_6 
                .shell_side_fluid_array_temperature()
                .unwrap();

            let ihx_shell_side_temp_profile_degc: Vec<f64> = 
                ihx_shell_side_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            let ihx_tube_side_temp_profile: Vec<ThermodynamicTemperature> = 
                ihx_sthe_6 
                .inner_tube_fluid_array_temperature()
                .unwrap();

            let ihx_tube_side_temp_profile_degc: Vec<f64> = 
                ihx_tube_side_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // steam generator tube side temp profile
            let sg_shell_side_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_steam_generator_shell_side_14 
                .pipe_fluid_array_temperature()
                .unwrap();

            let sg_shell_side_temp_profile_degc: Vec<f64> = 
                sg_shell_side_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // pipe 4, after reactor outlet 
            let pipe_4_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_4 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_4_temp_profile_degc: Vec<f64> = 
                pipe_4_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();
            // pipe 5, just before STHE
            let pipe_5_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_5 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_5_temp_profile_degc: Vec<f64> = 
                pipe_5_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();
            // pipe 7, just after STHE
            let pipe_7_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_7 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_7_temp_profile_degc: Vec<f64> = 
                pipe_7_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();
            // pipe 8, just before pump
            let pipe_8_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_8 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_8_temp_profile_degc: Vec<f64> = 
                pipe_8_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // pipe 10, just after pump
            let pipe_10_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_10 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_10_temp_profile_degc: Vec<f64> = 
                pipe_10_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // pipe 11, just before reactor inlet
            let pipe_11_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_11 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_11_temp_profile_degc: Vec<f64> = 
                pipe_11_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();


            // pipe 12, just before STHE tube side
            let pipe_12_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_12 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_12_temp_profile_degc: Vec<f64> = 
                pipe_12_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // pipe 13, just before steam generator shell side
            let pipe_13_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_13 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_13_temp_profile_degc: Vec<f64> = 
                pipe_13_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // pipe 15, just after steam generator shell side
            let pipe_15_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_15 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_15_temp_profile_degc: Vec<f64> = 
                pipe_15_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();


            // pipe 17, just after steam generator shell side
            let pipe_17_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pipe_17 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pipe_17_temp_profile_degc: Vec<f64> = 
                pipe_17_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // pri pump
            let pump_9_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_pri_loop_pump_9 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pump_9_temp_profile_degc: Vec<f64> = 
                pump_9_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // intrmd pump
            let pump_16_temp_profile: Vec<ThermodynamicTemperature> = 
                fhr_intrmd_loop_pump_16 
                .pipe_fluid_array_temperature()
                .unwrap();

            let pump_16_temp_profile_degc: Vec<f64> = 
                pump_16_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            // downcomer_temp profile
            let downcomer_2_temp_profile: Vec<ThermodynamicTemperature> = 
                downcomer_pipe_2 
                .pipe_fluid_array_temperature()
                .unwrap();

            let downcomer_2_temp_profile_degc: Vec<f64> = 
                downcomer_2_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();

            let downcomer_3_temp_profile: Vec<ThermodynamicTemperature> = 
                downcomer_pipe_3 
                .pipe_fluid_array_temperature()
                .unwrap();

            let downcomer_3_temp_profile_degc: Vec<f64> = 
                downcomer_3_temp_profile
                .into_iter()
                .map(|temperature|{
                    (temperature.get::<degree_celsius>()*100.0).round()/100.0
                })
            .collect();


            let fhr_state = FHRThermalHydraulicsState {
                reactor_branch_flow,
                downcomer_branch_1_flow,
                downcomer_branch_2_flow,
                intermediate_heat_exchanger_branch_flow: pri_loop_intermediate_heat_exchanger_branch_flow,
                intrmd_loop_ihx_br_flow,
                intrmd_loop_steam_gen_br_flow,
                simulation_time,
                reactor_temp_profile_degc,
                ihx_shell_side_temp_profile_degc,
                ihx_tube_side_temp_profile_degc,
                sg_shell_side_temp_profile_degc,
                pipe_4_temp_profile_degc,
                pipe_5_temp_profile_degc,
                pipe_7_temp_profile_degc,
                pipe_8_temp_profile_degc,
                pump_9_temp_profile_degc,
                pipe_10_temp_profile_degc,
                pipe_11_temp_profile_degc,
                pipe_12_temp_profile_degc,
                pipe_13_temp_profile_degc,
                pipe_15_temp_profile_degc,
                pump_16_temp_profile_degc,
                pipe_17_temp_profile_degc,
                downcomer_2_temp_profile_degc,
                downcomer_3_temp_profile_degc,
            };

            // if one wants to monitor flow through the loop
            let debugging = false;
            if debugging {
                dbg!(&fhr_state);
            }
            return fhr_state;
        }



    pub fn calculate_thermal_hydraulics_loop(
        fhr_state: Arc<Mutex<FHRState>>){

        let thermal_hydraulics_timestep = Time::new::<second>(0.1);

        let fhr_state_clone = fhr_state.clone();
        // now, time controls 
        let loop_time = SystemTime::now();
        let mut current_simulation_time = Time::ZERO;

        // create components first
        let initial_temperature = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.lock().unwrap().core_outlet_temp_degc
        );
        let mut reactor_pipe_1 = new_reactor_vessel_pipe_1(initial_temperature);
        let mut downcomer_pipe_2 = new_downcomer_pipe_2(initial_temperature);
        let mut downcomer_pipe_3 = new_downcomer_pipe_3(initial_temperature);

        // pri loop branch (positive is in this order of flow)
        let mut fhr_pipe_11 = new_fhr_pipe_11(initial_temperature);
        let mut fhr_pipe_10 = new_fhr_pipe_10(initial_temperature);
        let mut fhr_pri_loop_pump_9 = new_fhr_pri_loop_pump_9(initial_temperature);
        let mut fhr_pipe_8 = new_fhr_pipe_8(initial_temperature);
        let mut fhr_pipe_7 = new_fhr_pipe_7(initial_temperature);
        // note that for HITEC, the temperature range is from 
        // 440-800K 
        // this is 167-527C
        // so intial temperature of 500C is ok
        let mut ihx_sthe_6 = new_ihx_sthe_6_version_1(initial_temperature);
        let mut fhr_pipe_5 = new_fhr_pipe_5(initial_temperature);
        let mut fhr_pipe_4 = new_fhr_pipe_4_ver_2(initial_temperature);


        let initial_temperature_intrmd_loop = 
            initial_temperature;
        // intermediate loop ihx side 
        // (excluding sthe)
        let mut fhr_pipe_17 = new_fhr_pipe_17(initial_temperature_intrmd_loop);
        let mut fhr_pipe_12 = new_fhr_pipe_12(initial_temperature_intrmd_loop);

        // intermediate loop steam generator side 
        let mut fhr_intrmd_loop_pump_16 = new_fhr_intermediate_loop_pump_16(
            initial_temperature_intrmd_loop);
        let mut fhr_pipe_15 = new_fhr_pipe_15(initial_temperature_intrmd_loop);
        let mut fhr_steam_generator_shell_side_14 
            = new_fhr_intermediate_loop_steam_generator_shell_side_14(
                initial_temperature_intrmd_loop);
        let mut fhr_pipe_13 = new_fhr_pipe_13(initial_temperature_intrmd_loop);


        // probably want to use fhr state
        let mut pri_loop_pump_pressure = Pressure::new::<kilopascal>(-10.0);
        let mut intrmd_loop_pump_pressure = Pressure::new::<kilopascal>(-10.0);

        // mixing nodes for pri loop 
        let mut bottom_mixing_node_pri_loop = 
            gfhr_bottom_mixing_node_pri_loop(initial_temperature);
        let mut top_mixing_node_pri_loop = 
            gfhr_top_mixing_node_pri_loop(initial_temperature);
        // mixing nodes for intermediate loop 
        let mut bottom_mixing_node_intrmd_loop = 
            gfhr_bottom_mixing_node_intrmd_loop(initial_temperature_intrmd_loop);
        let mut top_mixing_node_intrmd_loop = 
            gfhr_top_mixing_node_intrmd_loop(initial_temperature_intrmd_loop);

        // create initial mass flowrates 

        // start with some initial flow rates
        let (mut reactor_branch_flow, mut downcomer_branch_1_flow, 
            mut downcomer_branch_2_flow, mut intermediate_heat_exchanger_branch_flow,
            mut intrmd_loop_ihx_br_flow,
            mut intrmd_loop_steam_gen_br_flow)
            = four_branch_pri_and_intermediate_loop_fluid_mechanics_only(
                pri_loop_pump_pressure, 
                intrmd_loop_pump_pressure, 
                &reactor_pipe_1, 
                &downcomer_pipe_2, 
                &downcomer_pipe_3, 
                &fhr_pipe_11, 
                &fhr_pipe_10, 
                &fhr_pri_loop_pump_9, 
                &fhr_pipe_8, 
                &fhr_pipe_7, 
                &ihx_sthe_6, 
                &fhr_pipe_5, 
                &fhr_pipe_4, 
                &fhr_pipe_17, 
                &fhr_pipe_12, 
                &fhr_intrmd_loop_pump_16, 
                &fhr_pipe_15, 
                &fhr_steam_generator_shell_side_14, 
                &fhr_pipe_13,
            );

        // steam generator settings 
        let steam_generator_tube_side_temperature = 
            ThermodynamicTemperature::new::<degree_celsius>(30.0);

        // I made this based on UA for 35 MWth heat load, and 
        // 30 degrees steam temperature, 300 degrees salt temperature
        let steam_generator_overall_ua: ThermalConductance 
            = ThermalConductance::new::<watt_per_kelvin>(1.2e5);

        let mut current_fhr_thermal_hydraulics_state = FHRThermalHydraulicsState {
            downcomer_branch_1_flow,
            downcomer_branch_2_flow,
            intermediate_heat_exchanger_branch_flow,
            intrmd_loop_ihx_br_flow,
            intrmd_loop_steam_gen_br_flow,
            reactor_branch_flow,
            simulation_time: current_simulation_time,
            reactor_temp_profile_degc: vec![],
            ihx_shell_side_temp_profile_degc: vec![],
            ihx_tube_side_temp_profile_degc: vec![],
            sg_shell_side_temp_profile_degc: vec![],
            pipe_4_temp_profile_degc: vec![],
            pipe_5_temp_profile_degc: vec![],
            pipe_7_temp_profile_degc: vec![],
            pipe_8_temp_profile_degc: vec![],
            pump_9_temp_profile_degc: vec![],
            pipe_10_temp_profile_degc: vec![],
            pipe_11_temp_profile_degc: vec![],
            pipe_12_temp_profile_degc: vec![],
            pipe_13_temp_profile_degc: vec![],
            pipe_15_temp_profile_degc: vec![],
            pump_16_temp_profile_degc: vec![],
            pipe_17_temp_profile_degc: vec![],
            downcomer_2_temp_profile_degc: vec![],
            downcomer_3_temp_profile_degc: vec![],
        };
        dbg!(&current_fhr_thermal_hydraulics_state);
        // calculation loop (indefinite)
        //
        // to be done once every timestep
        loop {

            // so now, let's do the necessary things
            // first, timestep and loop time 
            //
            // second, read and update the local_ciet_state

            let loop_time_start = loop_time.elapsed().unwrap();

            //
            let accumulated_energy_from_prke = 
                Energy::new::<kilojoule>(
                    fhr_state_clone.lock().unwrap().prke_loop_accumulated_heat_removal_kilojoules
                );

            // once energy is "taken out" from the store, 
            // it is meant to be zero
            fhr_state_clone.lock().unwrap().prke_loop_accumulated_heat_removal_kilojoules 
                = 0.0;
            let reactor_power = 
                accumulated_energy_from_prke/thermal_hydraulics_timestep;

            // then primary loop pressure 
            
            pri_loop_pump_pressure = 
                Pressure::new::<kilopascal>(
                    -fhr_state_clone.lock().unwrap().fhr_pri_loop_pump_pressure_kilopascals
                );
            intrmd_loop_pump_pressure = 
                Pressure::new::<kilopascal>(
                    -fhr_state_clone.lock().unwrap().fhr_intermediate_loop_pump_pressure_kilopascals
                );

            current_fhr_thermal_hydraulics_state = 
                Self::four_branch_pri_and_intermediate_loop_single_time_step(
                    pri_loop_pump_pressure, 
                    intrmd_loop_pump_pressure, 
                    reactor_power, 
                    thermal_hydraulics_timestep, current_simulation_time, 
                    &mut reactor_pipe_1, &mut downcomer_pipe_2, 
                    &mut downcomer_pipe_3, &mut bottom_mixing_node_pri_loop, 
                    &mut top_mixing_node_pri_loop, &mut fhr_pipe_11, 
                    &mut fhr_pipe_10, &mut fhr_pri_loop_pump_9, 
                    &mut fhr_pipe_8, &mut fhr_pipe_7, &mut ihx_sthe_6, 
                    &mut fhr_pipe_5, &mut fhr_pipe_4, &mut fhr_pipe_17, 
                    &mut fhr_pipe_12, &mut fhr_intrmd_loop_pump_16, 
                    &mut fhr_pipe_15, &mut fhr_steam_generator_shell_side_14, 
                    &mut fhr_pipe_13, &mut bottom_mixing_node_intrmd_loop, 
                    &mut top_mixing_node_intrmd_loop, 
                    steam_generator_tube_side_temperature, 
                    steam_generator_overall_ua);

            dbg!(&current_fhr_thermal_hydraulics_state);

            current_simulation_time += thermal_hydraulics_timestep;

            let simulation_time_seconds = current_simulation_time.get::<second>();

            let elapsed_time_seconds = 
                (loop_time.elapsed().unwrap().as_secs_f64() * 100.0).round()/100.0;

            *&mut fhr_state_clone.lock().unwrap().thermal_hydraulics_simulation_time_seconds 
                = elapsed_time_seconds;

            let overall_simulation_in_realtime_or_faster: bool = 
                simulation_time_seconds > elapsed_time_seconds;

            // now update the ciet state 
            let loop_time_end = loop_time.elapsed().unwrap();
            let time_taken_for_calculation_loop_microseconds: f64 = 
                (loop_time_end - loop_time_start)
                .as_micros() as f64;

            *&mut fhr_state_clone.lock().unwrap().thermal_hydraulics_timestep_microseconds
                = thermal_hydraulics_timestep.get::<microsecond>().round();
            *&mut fhr_state_clone.lock().unwrap().thermal_hydraulics_calc_time_microseconds
                = time_taken_for_calculation_loop_microseconds;

            // update temperatures
            {
                let mut fhr_state_lock = fhr_state_clone.lock().unwrap();

                // the reactor branch itself has five elements
                // from 0,1,2,3,4 in order of going from bottom 
                // of the core to the top
                fhr_state_lock.core_inlet_temp_degc = 
                    current_fhr_thermal_hydraulics_state
                    .reactor_temp_profile_degc[0];
                fhr_state_lock.core_bottom_temp_degc = 
                    current_fhr_thermal_hydraulics_state
                    .reactor_temp_profile_degc[1];
                fhr_state_lock.pebble_bed_coolant_temp_degc = 
                    current_fhr_thermal_hydraulics_state
                    .reactor_temp_profile_degc[2];
                fhr_state_lock.core_top_temp_degc = 
                    current_fhr_thermal_hydraulics_state
                    .reactor_temp_profile_degc[3];
                fhr_state_lock.core_outlet_temp_degc = 
                    current_fhr_thermal_hydraulics_state
                    .reactor_temp_profile_degc[4];

                // the downcomers 1 and 2 also have branches
                // with 5 nodes each 
                // however, the fhr itself only has three distinct regions 
                // for display
                // for now I will just use node 0, 2 and 4
                //

                fhr_state_lock.left_downcomer_lower_temp_degc = 
                    current_fhr_thermal_hydraulics_state 
                    .downcomer_2_temp_profile_degc[0];
                fhr_state_lock.left_downcomer_mid_temp_degc = 
                    current_fhr_thermal_hydraulics_state 
                    .downcomer_2_temp_profile_degc[2];
                fhr_state_lock.left_downcomer_upper_temp_degc = 
                    current_fhr_thermal_hydraulics_state 
                    .downcomer_2_temp_profile_degc[4];

                fhr_state_lock.right_downcomer_lower_temp_degc = 
                    current_fhr_thermal_hydraulics_state 
                    .downcomer_3_temp_profile_degc[0];
                fhr_state_lock.right_downcomer_mid_temp_degc = 
                    current_fhr_thermal_hydraulics_state 
                    .downcomer_3_temp_profile_degc[2];
                fhr_state_lock.right_downcomer_upper_temp_degc = 
                    current_fhr_thermal_hydraulics_state 
                    .downcomer_3_temp_profile_degc[4];

                // flowrate diagnostics 
                fhr_state_lock.intermediate_loop_clockwise_flow_kg_per_s = 
                    (current_fhr_thermal_hydraulics_state 
                    .intrmd_loop_ihx_br_flow 
                    .get::<kilogram_per_second>()*1000.0)/1000.0;

                fhr_state_lock.reactor_branch_flowrate_kg_per_s = 
                    (current_fhr_thermal_hydraulics_state 
                    .reactor_branch_flow 
                    .get::<kilogram_per_second>()*1000.0)/1000.0;
                fhr_state_lock.downcomer1_branch_flowrate_kg_per_s = 
                    (current_fhr_thermal_hydraulics_state 
                    .downcomer_branch_1_flow 
                    .get::<kilogram_per_second>()*1000.0)/1000.0;

                fhr_state_lock.downcomer2_branch_flowrate_kg_per_s = 
                    (current_fhr_thermal_hydraulics_state 
                    .downcomer_branch_2_flow
                    .get::<kilogram_per_second>()*1000.0)/1000.0;
                fhr_state_lock.ihx_branch_flowrate_kg_per_s = 
                    (current_fhr_thermal_hydraulics_state 
                    .intermediate_heat_exchanger_branch_flow
                    .get::<kilogram_per_second>()*1000.0)/1000.0;

            }


            let time_to_sleep_microseconds: u64 = 
                (thermal_hydraulics_timestep.get::<microsecond>() - 
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

            


        }

    }
}

/// contains simple components for the fhr simulator
///
/// these are components for primary loop and secondary loop 
/// turbine components not included (will be in tampines-steam-tables)
pub mod components;

/// contains functions for calculating pri loop 
/// fluid mechanics
pub mod pri_loop_fluid_mechanics_calc_fns;


/// code for fhr thermal hydraulics state 
pub mod fhr_thermal_hydraulics_state;
