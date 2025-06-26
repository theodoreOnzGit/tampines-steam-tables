/// this is the struct to contain the thermal hydraulics state of the fhr
use uom::si::f64::*;
#[derive(Debug,Clone)]
pub(crate) struct FHRThermalHydraulicsState {
    /// reactor branch flow (upwards through the core)
    /// note that positive flow means from bottom mixing node to top
    pub reactor_branch_flow: MassRate,
    /// downcomer 1 branch flow (upwards through the core)
    /// note that positive flow means from bottom mixing node to top
    pub downcomer_branch_1_flow: MassRate,
    /// downcomer 2 branch flow (upwards through the core)
    /// note that positive flow means from bottom mixing node to top
    pub downcomer_branch_2_flow: MassRate,
    /// ihx branch flow 
    /// note that positive flow means from bottom mixing node to top
    pub intermediate_heat_exchanger_branch_flow: MassRate,
    /// ihx branch flow 
    /// note that positive flow means from bottom 
    /// (between pipe 17 and pump 16) 
    /// to top
    /// (between pipe 12 and pipe 13)
    pub intrmd_loop_ihx_br_flow: MassRate,
    /// steam generator branch
    /// note that positive flow means from bottom 
    /// (between pipe 17 and pump 16) 
    /// to top
    /// (between pipe 12 and pipe 13)
    pub intrmd_loop_steam_gen_br_flow: MassRate,

    // other diagnostics 
    /// shows the current simulation time
    pub simulation_time: Time,

    // temperature diagnostics 
    /// shows the current reactor temperature profile in degc (2dp)
    pub reactor_temp_profile_degc: Vec<f64>,
    /// shows the current ihx shell side temperature profile in degc (2dp)
    pub ihx_shell_side_temp_profile_degc: Vec<f64>,
    /// shows the current ihx tube side temperature profile in degc (2dp)
    pub ihx_tube_side_temp_profile_degc: Vec<f64>,
    /// shows the current steam generator side temperature profile in degc (2dp)
    pub sg_shell_side_temp_profile_degc: Vec<f64>,

    /// shows the temperature profile of pipe_4
    pub pipe_4_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_5
    pub pipe_5_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_7
    pub pipe_7_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_8
    pub pipe_8_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pump_9 in the primary loop
    pub pump_9_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_10
    pub pipe_10_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_11
    pub pipe_11_temp_profile_degc: Vec<f64>,


    // intermediate loop

    /// shows the temperature profile of pipe_12
    pub pipe_12_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_13
    pub pipe_13_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_15
    pub pipe_15_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pump_16 in the intermediate loop
    pub pump_16_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_17
    pub pipe_17_temp_profile_degc: Vec<f64>,

    // downcomers
    /// shows the temperature profile of pipe_12
    pub downcomer_2_temp_profile_degc: Vec<f64>,
    /// shows the temperature profile of pipe_13
    pub downcomer_3_temp_profile_degc: Vec<f64>,

}

