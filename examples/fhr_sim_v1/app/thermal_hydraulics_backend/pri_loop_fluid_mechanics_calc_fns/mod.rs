use std::ops::Deref;
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;

use tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component::FluidComponent;
use tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollection;
use tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollectionMethods;
use tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_super_collection::FluidComponentSuperCollection;
use tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidComponentTrait;
use tuas_boussinesq_solver::pre_built_components::ciet_three_branch_plus_dracs::solver_functions::get_mass_flowrate_two_branches;
use tuas_boussinesq_solver::pre_built_components::ciet_three_branch_plus_dracs::solver_functions::get_mass_flowrate_vector_for_dhx_heater_and_ctah_branches;
use tuas_boussinesq_solver::pre_built_components::insulated_pipes_and_fluid_components::InsulatedFluidComponent;
use tuas_boussinesq_solver::pre_built_components::non_insulated_fluid_components::NonInsulatedFluidComponent;
use tuas_boussinesq_solver::pre_built_components::shell_and_tube_heat_exchanger::SimpleShellAndTubeHeatExchanger;
use tuas_boussinesq_solver::prelude::beta_testing::InsulatedPorousMediaFluidComponent;
use uom::ConstZero;
use uom::si::f64::*;


pub mod parallel_branch_flow_calculator;

/// for the gFHR primary loop, and intermediate loop 
/// there are four branches that need to be solved for flowrate 
///
/// this code handles the solution procedure
/// using the tuas_boussinesq_solver library code
///
/// and only handles fluid mechanics (isothermal)
pub fn four_branch_pri_and_intermediate_loop_fluid_mechanics_only(
    pri_loop_pump_pressure: Pressure,
    intrmd_loop_pump_pressure: Pressure,
    // reactor branch
    reactor_pipe_1: &InsulatedFluidComponent,
    // downcomer branch 1
    downcomer_pipe_2: &InsulatedFluidComponent,
    // downcomer branch 2
    downcomer_pipe_3: &InsulatedFluidComponent,
    // Intermediate heat exchanger branch in pri loop
    fhr_pipe_11: &InsulatedFluidComponent,
    fhr_pipe_10: &InsulatedFluidComponent,
    fhr_pri_loop_pump_9: &NonInsulatedFluidComponent,
    fhr_pipe_8: &InsulatedFluidComponent,
    fhr_pipe_7: &InsulatedFluidComponent,
    ihx_sthe_6: &SimpleShellAndTubeHeatExchanger,
    fhr_pipe_5: &InsulatedFluidComponent,
    fhr_pipe_4: &InsulatedFluidComponent,
    // intermediate loop ihx side
    fhr_pipe_17: &InsulatedFluidComponent,
    fhr_pipe_12: &InsulatedFluidComponent,
    // intermediate loop steam generator side
    fhr_intrmd_loop_pump_16: &NonInsulatedFluidComponent,
    fhr_pipe_15: &InsulatedFluidComponent,
    fhr_steam_generator_shell_side_14: &NonInsulatedFluidComponent,
    fhr_pipe_13: &InsulatedFluidComponent,

    ) -> (MassRate, MassRate, MassRate, MassRate, MassRate, MassRate){

        // pri loop

        let mut reactor_branch = 
            FluidComponentCollection::new_series_component_collection();

        reactor_branch.clone_and_add_component(reactor_pipe_1);




        let mut pri_downcomer_branch_1 = 
            FluidComponentCollection::new_series_component_collection();

        pri_downcomer_branch_1.clone_and_add_component(downcomer_pipe_2);




        let mut pri_downcomer_branch_2 = 
            FluidComponentCollection::new_series_component_collection();

        pri_downcomer_branch_2.clone_and_add_component(downcomer_pipe_3);




        let mut pri_loop_intermediate_heat_exchanger_branch =
            FluidComponentCollection::new_series_component_collection();

        let mut fhr_pri_loop_pump_9_with_pressure_set = fhr_pri_loop_pump_9.clone();
        fhr_pri_loop_pump_9_with_pressure_set.set_internal_pressure_source(pri_loop_pump_pressure);
        let ihx_shell_side_6_clone = ihx_sthe_6.get_clone_of_shell_side_fluid_component();

        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_11);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_10);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            &fhr_pri_loop_pump_9_with_pressure_set);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_8);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_7);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            &ihx_shell_side_6_clone);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_5);
        pri_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_4);




        // intermediate loop
        // ihx side

        let mut intrmd_loop_intermediate_heat_exchanger_branch =
            FluidComponentCollection::new_series_component_collection();

        intrmd_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_17);
        intrmd_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            &ihx_sthe_6.get_clone_of_tube_side_parallel_tube_fluid_component());
        intrmd_loop_intermediate_heat_exchanger_branch.clone_and_add_component(
            fhr_pipe_12);

        // intermediate loop
        // steam generator side

        let mut intrmd_loop_steam_generator_branch =
            FluidComponentCollection::new_series_component_collection();

        let mut fhr_intrmd_loop_pump_16_with_pressure_set = 
            fhr_intrmd_loop_pump_16.clone();
        fhr_intrmd_loop_pump_16_with_pressure_set
            .set_internal_pressure_source(intrmd_loop_pump_pressure);
        intrmd_loop_steam_generator_branch.clone_and_add_component(
            &fhr_intrmd_loop_pump_16_with_pressure_set);
        intrmd_loop_steam_generator_branch.clone_and_add_component(
            fhr_pipe_15);
        intrmd_loop_steam_generator_branch.clone_and_add_component(
            fhr_steam_generator_shell_side_14);
        intrmd_loop_steam_generator_branch.clone_and_add_component(
            fhr_pipe_13);


        // calculate pri loop side fluid mechanics
        let mut pri_loop_branches = 
            FluidComponentSuperCollection::default();

        pri_loop_branches.set_orientation_to_parallel();

        pri_loop_branches.fluid_component_super_vector.push(reactor_branch);
        pri_loop_branches.fluid_component_super_vector.push(pri_downcomer_branch_1);
        pri_loop_branches.fluid_component_super_vector.push(pri_downcomer_branch_2);

        pri_loop_branches.fluid_component_super_vector.push(pri_loop_intermediate_heat_exchanger_branch);

        let pressure_change_across_each_branch_pri_loop = 
            pri_loop_branches.get_pressure_change(MassRate::ZERO);


        let pri_loop_mass_rate_vector 
            = pri_loop_branches.get_mass_flowrate_across_each_parallel_branch(
                pressure_change_across_each_branch_pri_loop);

        let (reactor_branch_flow, downcomer_branch_1_flow,
            downcomer_branch_2_flow, intermediate_heat_exchanger_branch_flow)
            = (pri_loop_mass_rate_vector[0], pri_loop_mass_rate_vector[1],
                pri_loop_mass_rate_vector[2], pri_loop_mass_rate_vector[3]);

        // calculate intermediate loop side fluid mechanics

        let mut intrmd_loop_branches = 
            FluidComponentSuperCollection::default();

        intrmd_loop_branches.set_orientation_to_parallel();

        intrmd_loop_branches.fluid_component_super_vector.push(
            intrmd_loop_intermediate_heat_exchanger_branch);
        intrmd_loop_branches.fluid_component_super_vector.push(
            intrmd_loop_steam_generator_branch);

        let pressure_change_across_each_branch_intrmd_loop = 
            intrmd_loop_branches.get_pressure_change(MassRate::ZERO);

        let intrmd_loop_mass_rate_vector 
            = intrmd_loop_branches.get_mass_flowrate_across_each_parallel_branch(
                pressure_change_across_each_branch_intrmd_loop);
        let (intrmd_loop_ihx_br_flow, intrmd_loop_steam_gen_br_flow) = 
            (intrmd_loop_mass_rate_vector[0],
             intrmd_loop_mass_rate_vector[1]);


        return (reactor_branch_flow, downcomer_branch_1_flow,
            downcomer_branch_2_flow, intermediate_heat_exchanger_branch_flow,
            intrmd_loop_ihx_br_flow, intrmd_loop_steam_gen_br_flow);
}

/// fluid mechanics bit for primary loop 
/// calculate fluid mechanics across four branches in parallel,
/// ie  the mass flowrate across each branch
///
/// this is for the fhr having 
/// 1. reactor branch
/// 2. downcomer branch 1 
/// 3. downcomer branch 2 
/// 4. intermediate heat exchanger branch
pub fn get_mass_flowrate_across_for_reactor_downcomers_and_heat_exchg_br(
    pri_loop_branches: &FluidComponentSuperCollection
) -> (MassRate, MassRate, MassRate, MassRate) {

    // basically the net flowrate across the branches as a whole is zero 
    //
    let pressure_change_across_each_branch = 
        pri_loop_branches.get_pressure_change(MassRate::ZERO);

    let mass_flowrate_across_each_branch: Vec<MassRate> = 
        pri_loop_branches. 
        get_mass_flowrate_across_each_parallel_branch(
            pressure_change_across_each_branch
        );

    // note, the mass flowrate order depends on how you add the branches 

    let mass_flowrate_branch_1 = mass_flowrate_across_each_branch[0];
    let mass_flowrate_branch_2 = mass_flowrate_across_each_branch[1];
    let mass_flowrate_branch_3 = mass_flowrate_across_each_branch[2];
    let mass_flowrate_branch_4 = mass_flowrate_across_each_branch[3];


    return(mass_flowrate_branch_1,
        mass_flowrate_branch_2,
        mass_flowrate_branch_3,
        mass_flowrate_branch_4);


}


/// ciet coupled dracs loop calculations, serves as an example
pub fn ciet_get_abs_mass_flowrate_across_two_branches(dracs_branches: &FluidComponentSuperCollection) -> 
MassRate {
    let pressure_change_across_each_branch = 
        dracs_branches.get_pressure_change(MassRate::ZERO);

    let mass_flowrate_across_each_branch: Vec<MassRate> = 
        dracs_branches.
        get_mass_flowrate_across_each_parallel_branch(
            pressure_change_across_each_branch
        );

    let mut mass_flowrate: MassRate = 
        mass_flowrate_across_each_branch[0];


    // get absolute value
    mass_flowrate = mass_flowrate.abs();

    mass_flowrate

}
/// ciet coupled dracs loop calculations, serves as an example
pub fn ciet_coupled_dracs_fluid_mechanics_calc_abs_mass_rate_no_tchx_calibration(
    pipe_34: &InsulatedFluidComponent,
    pipe_33: &InsulatedFluidComponent,
    pipe_32: &InsulatedFluidComponent,
    pipe_31a: &InsulatedFluidComponent,
    static_mixer_61_label_31: &InsulatedFluidComponent,
    dhx_tube_side_30b: &NonInsulatedFluidComponent,
    dhx_tube_side_heat_exchanger_30: &FluidComponent,
    dhx_tube_side_30a: &NonInsulatedFluidComponent,
    tchx_35a: &NonInsulatedFluidComponent,
    tchx_35b: &NonInsulatedFluidComponent,
    static_mixer_60_label_36: &InsulatedFluidComponent,
    pipe_36a: &InsulatedFluidComponent,
    pipe_37: &InsulatedFluidComponent,
    flowmeter_60_37a: &NonInsulatedFluidComponent,
    pipe_38: &InsulatedFluidComponent,
    pipe_39: &InsulatedFluidComponent,
)-> MassRate {

    let mut dracs_hot_branch = 
        FluidComponentCollection::new_series_component_collection();

    dracs_hot_branch.clone_and_add_component(pipe_34);
    dracs_hot_branch.clone_and_add_component(pipe_33);
    dracs_hot_branch.clone_and_add_component(pipe_32);
    dracs_hot_branch.clone_and_add_component(pipe_31a);
    dracs_hot_branch.clone_and_add_component(static_mixer_61_label_31);
    dracs_hot_branch.clone_and_add_component(dhx_tube_side_30b);
    dracs_hot_branch.clone_and_add_component(dhx_tube_side_heat_exchanger_30);
    dracs_hot_branch.clone_and_add_component(dhx_tube_side_30a);


    let mut dracs_cold_branch = 
        FluidComponentCollection::new_series_component_collection();

    dracs_cold_branch.clone_and_add_component(tchx_35a);
    dracs_cold_branch.clone_and_add_component(tchx_35b);
    dracs_cold_branch.clone_and_add_component(static_mixer_60_label_36);
    dracs_cold_branch.clone_and_add_component(pipe_36a);
    dracs_cold_branch.clone_and_add_component(pipe_37);
    dracs_cold_branch.clone_and_add_component(flowmeter_60_37a);
    dracs_cold_branch.clone_and_add_component(pipe_38);
    dracs_cold_branch.clone_and_add_component(pipe_39);

    let mut dracs_branches = 
        FluidComponentSuperCollection::default();

    dracs_branches.set_orientation_to_parallel();
    dracs_branches.fluid_component_super_vector.push(dracs_hot_branch);
    dracs_branches.fluid_component_super_vector.push(dracs_cold_branch);

    let abs_mass_rate = ciet_get_abs_mass_flowrate_across_two_branches(&dracs_branches);

    abs_mass_rate

}


pub fn ciet_three_branch_pri_loop_flowrates_parallel_ver_4(
    pump_pressure: Pressure,
    ctah_branch_blocked: bool,
    dhx_branch_blocked: bool,
    pipe_4: &InsulatedFluidComponent,
    pipe_3: &InsulatedFluidComponent,
    pipe_2a: &InsulatedFluidComponent,
    static_mixer_10_label_2: &InsulatedFluidComponent,
    heater_top_head_1a: &InsulatedPorousMediaFluidComponent,
    heater_ver_1: &InsulatedPorousMediaFluidComponent,
    heater_bottom_head_1b: &InsulatedPorousMediaFluidComponent,
    pipe_18: &InsulatedFluidComponent,
    pipe_5a: &InsulatedFluidComponent,
    pipe_26: &InsulatedFluidComponent,
    pipe_25a: &InsulatedFluidComponent,
    static_mixer_21_label_25: &InsulatedFluidComponent,
    dhx_shell_side_pipe_24: &FluidComponent,
    static_mixer_20_label_23: &InsulatedFluidComponent,
    pipe_23a: &InsulatedFluidComponent,
    pipe_22: &InsulatedFluidComponent,
    flowmeter_20_21a: &NonInsulatedFluidComponent,
    pipe_21: &InsulatedFluidComponent,
    pipe_20: &InsulatedFluidComponent,
    pipe_19: &InsulatedFluidComponent,
    pipe_17b: &InsulatedFluidComponent,
    // ctah branch
    pipe_5b: &InsulatedFluidComponent,
    static_mixer_41_label_6 :&InsulatedFluidComponent,
    pipe_6a :&InsulatedFluidComponent,
    ctah_vertical_label_7a :&NonInsulatedFluidComponent,
    ctah_horizontal_label_7b :&NonInsulatedFluidComponent,
    pipe_8a :&InsulatedFluidComponent,
    static_mixer_40_label_8 :&InsulatedFluidComponent,
    pipe_9 :&InsulatedFluidComponent,
    pipe_10 :&InsulatedFluidComponent,
    pipe_11 :&InsulatedFluidComponent,
    pipe_12 :&InsulatedFluidComponent,
    ctah_pump :&NonInsulatedFluidComponent,
    pipe_13 : &InsulatedFluidComponent,
    pipe_14 : &InsulatedFluidComponent,
    flowmeter_40_14a :&NonInsulatedFluidComponent,
    pipe_15 :&InsulatedFluidComponent,
    pipe_16 :&InsulatedFluidComponent,
    pipe_17a :&InsulatedFluidComponent,
    ) ->
(MassRate, MassRate, MassRate) {


    let mut heater_branch = 
        FluidComponentCollection::new_series_component_collection();

    heater_branch.clone_and_add_component(pipe_4);
    heater_branch.clone_and_add_component(pipe_3);
    heater_branch.clone_and_add_component(pipe_2a);
    heater_branch.clone_and_add_component(static_mixer_10_label_2);
    heater_branch.clone_and_add_component(heater_top_head_1a);
    heater_branch.clone_and_add_component(heater_ver_1);
    heater_branch.clone_and_add_component(heater_bottom_head_1b);
    heater_branch.clone_and_add_component(pipe_18);


    let mut dhx_branch = 
        FluidComponentCollection::new_series_component_collection();

    dhx_branch.clone_and_add_component(pipe_5a);
    dhx_branch.clone_and_add_component(pipe_26);
    dhx_branch.clone_and_add_component(pipe_25a);
    dhx_branch.clone_and_add_component(static_mixer_21_label_25);
    dhx_branch.clone_and_add_component(dhx_shell_side_pipe_24);
    dhx_branch.clone_and_add_component(static_mixer_20_label_23);
    dhx_branch.clone_and_add_component(pipe_23a);
    dhx_branch.clone_and_add_component(pipe_22);
    dhx_branch.clone_and_add_component(flowmeter_20_21a);
    dhx_branch.clone_and_add_component(pipe_21);
    dhx_branch.clone_and_add_component(pipe_20);
    dhx_branch.clone_and_add_component(pipe_19);
    dhx_branch.clone_and_add_component(pipe_17b);

    let mut ctah_branch = FluidComponentCollection::new_series_component_collection();
    
    ctah_branch.clone_and_add_component(pipe_5b);
    ctah_branch.clone_and_add_component(static_mixer_41_label_6);
    ctah_branch.clone_and_add_component(pipe_6a);
    ctah_branch.clone_and_add_component(ctah_vertical_label_7a);
    ctah_branch.clone_and_add_component(ctah_horizontal_label_7b);
    ctah_branch.clone_and_add_component(pipe_8a);
    ctah_branch.clone_and_add_component(static_mixer_40_label_8);
    ctah_branch.clone_and_add_component(pipe_9);
    ctah_branch.clone_and_add_component(pipe_10);
    ctah_branch.clone_and_add_component(pipe_11);
    ctah_branch.clone_and_add_component(pipe_12);
    let mut ctah_pump_clone: NonInsulatedFluidComponent 
        = ctah_pump.clone();
    ctah_pump_clone.set_internal_pressure_source(pump_pressure);
    ctah_branch.clone_and_add_component(&ctah_pump_clone);
    ctah_branch.clone_and_add_component(pipe_13);
    ctah_branch.clone_and_add_component(pipe_14);
    ctah_branch.clone_and_add_component(flowmeter_40_14a);
    ctah_branch.clone_and_add_component(pipe_15);
    ctah_branch.clone_and_add_component(pipe_16);
    ctah_branch.clone_and_add_component(pipe_17a);

    let mut pri_loop_branches = 
        FluidComponentSuperCollection::default();

    pri_loop_branches.set_orientation_to_parallel();

    if ctah_branch_blocked && !dhx_branch_blocked {

        pri_loop_branches.fluid_component_super_vector.push(dhx_branch);
        pri_loop_branches.fluid_component_super_vector.push(heater_branch);

        let (dhx_flow, heater_flow) = 
            get_mass_flowrate_two_branches(
                &pri_loop_branches);

        return (dhx_flow, heater_flow, MassRate::ZERO);

    } else if dhx_branch_blocked && !ctah_branch_blocked {

        pri_loop_branches.fluid_component_super_vector.push(heater_branch);
        pri_loop_branches.fluid_component_super_vector.push(ctah_branch);

        let (heater_flow, ctah_flow) = 
            get_mass_flowrate_two_branches(
                &pri_loop_branches);
        return (MassRate::ZERO, heater_flow, ctah_flow);

    } else if ctah_branch_blocked && dhx_branch_blocked {
        // all flows blocked, no nothing to see here

        return (MassRate::ZERO, MassRate::ZERO, MassRate::ZERO);
    } else {
        
        // two loops open scenario 

        let mut pri_loop_ctah_and_heater_br =
            FluidComponentSuperCollection::default();

        pri_loop_ctah_and_heater_br.fluid_component_super_vector.push(heater_branch.clone());
        pri_loop_ctah_and_heater_br.fluid_component_super_vector.push(ctah_branch.clone());

        // create pointers for two branch flow 

        let heater_br_flow_two_br_scenario = Arc::new(Mutex::new(
                MassRate::ZERO));
        let ctah_br_flow_two_br_scenario = Arc::new(Mutex::new(
                MassRate::ZERO));

        let heater_br_flow_two_br_scenario_clone = heater_br_flow_two_br_scenario.clone();
        let ctah_br_flow_two_br_scenario_clone = ctah_br_flow_two_br_scenario.clone();

        // spawn a thread 

        let two_branch_flow_scenario_join_handle = 
            thread::spawn(move ||{
                let (heater_flow_two_br, ctah_flow_two_br) = 
                    get_mass_flowrate_two_branches(
                        &pri_loop_ctah_and_heater_br);

                // load this into the arc mutex ptr 

                *heater_br_flow_two_br_scenario_clone.lock().unwrap().deref_mut()
                    = heater_flow_two_br;
                *ctah_br_flow_two_br_scenario_clone.lock().unwrap().deref_mut()
                    = ctah_flow_two_br;
            });

        // all loops opened scenario
        pri_loop_branches.fluid_component_super_vector.push(dhx_branch);
        pri_loop_branches.fluid_component_super_vector.push(heater_branch);
        pri_loop_branches.fluid_component_super_vector.push(ctah_branch);

        let (dhx_flow, heater_flow, ctah_flow) = 
            get_mass_flowrate_vector_for_dhx_heater_and_ctah_branches(
                &pri_loop_branches);



        // if dhx flow is downwards, (positive flow, is ok)
        // if negative flow, then block it 

        let flow_diode_block_flow: bool = dhx_flow < MassRate::ZERO;

        // load the two branch solution
        two_branch_flow_scenario_join_handle.join().unwrap();

        if flow_diode_block_flow {


            let dhx_flow = MassRate::ZERO;
            let heater_flow: MassRate = *heater_br_flow_two_br_scenario.lock().unwrap().deref();
            let ctah_flow: MassRate = *ctah_br_flow_two_br_scenario.lock().unwrap().deref();
            return (dhx_flow, heater_flow, ctah_flow);

        } else {

            return (dhx_flow, heater_flow, ctah_flow);
        }


    }

}
