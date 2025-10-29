
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release
/// this represents the first iteration 
/// of the fhr simulator
///
/// basically one can do a FHR loop 
/// with a permenantly steady state steam cycle
/// the latter uses the tampines-steam-tables
fn main(){

    fhr_simulator_v1().unwrap();


}
use std::{sync::{Arc, Mutex}, thread};

use uom::si::{f64::*, power::kilowatt};

use crate::app::{graph_data::PagePlotData, panel_enum::Panel};

/// this is how the fhr simulator runs 
///
/// there are a number of things to improve though, some that make this 
/// simulator a little annoying to use 
///
/// 1. Salt freezes and temperature overrides. Sometimes the salt temperature 
/// goes too low for the loop. In that case, the salt freezes. 
/// But programmatically, the loop thermal hydraulics cycle crashes. 
/// I'd rather just have the flow stop and give the user an option to 
/// that any said salt loop. Thus resetting all the temperatures in 
/// the loop.
///
/// 2. Salt overheating in the HITEC loop. 
/// This is where the salt in HITEC loop gets too hot, thus making the 
/// HITEC thermally decompose in real life. Programmatically, the thermal 
/// hydraulics loop stops working because the temperature is out of range.
/// I'd rather give the user an option to restart the loop completely.
///
/// 3. For natural circulation, the loop flowrates sometimes don't add 
/// up to 0.0 kg/s. 
///
/// 4. For natural circulation, one should be able to block the downcomer 
/// valve in the primary loop. Thus having natural circulation 
///
/// 5. I'd rather have flowrates a little bit higher in the loop, getting to 
/// levels nearer the KP-FHR.
///
/// 6. The steam cycle is steady state yes. But this often causes numerical 
/// instabilities where the steam temperature is sometimes too high. 
/// This is especially when the UA value exceeds some number. It would be 
/// good to have some realism added where UA is not toggled by the user,
/// but is controlled via some departure from nucleate boiling model. In  
/// this manner, the steam shouldn't overheat unnaturally. 
///
/// In the simplest assumption, pool boiling can be assumed
/// https://www.nuclear-power.com/nuclear-engineering/heat-transfer/boiling-and-condensation/boiling-crisis-critical-heat-flux/
///
/// 7. Aesthetics. These are horrible for now and will need to change.
/// Long list to try and improve for later.
///
/// 8. Data presentation. 
///
/// 9. Validation and Verification. For PRKE, this is still not done,
/// will need to ensure that PRKE solver is reasonably accurate.
///
/// 
///
///
pub fn fhr_simulator_v1() -> eframe::Result<()> {
    env_logger::init(); // Log to stderr (if you run with `RUST_LOG=debug`).

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([800.0, 800.0]),
        ..Default::default()
    };
    eframe::run_native(
        "FHR Core / Primary Simulator V1 Powered by TUAS and teh-o-prke",
        native_options,
        Box::new(|cc| {
            // image support,
            // from 
            // https://github.com/emilk/egui/tree/master/examples/images
            egui_extras::install_image_loaders(&cc.egui_ctx);
            Ok(Box::new(FHRSimulatorApp::new(cc)))

    }

        ),
    )
}
#[derive(serde::Deserialize, serde::Serialize)]
#[serde(default)] // if we add new fields, give them default values when deserializing old state
#[derive(Clone, Debug)]
pub struct FHRSimulatorApp {

    pub fhr_state: Arc<Mutex<FHRState>>,

    /// what panel is open
    pub open_panel: Panel,

    #[serde(skip)]
    /// pointer for plotting 
    pub fhr_simulator_ptr_for_plotting: Arc<Mutex<PagePlotData>>
}

#[derive(serde::Deserialize, serde::Serialize)]
#[serde(default)] // if we add new fields, give them default values when deserializing old state
#[derive(Clone,Debug)]
pub struct FHRState {
    /// left control rod insertion fraction
    pub left_cr_insertion_frac: f32,
    /// right control rod insertion fraction
    pub right_cr_insertion_frac: f32,

    // temperatures for both reactor feedback and display
    pub pebble_core_temp_degc: f64,
    pub pebble_bed_coolant_temp_degc: f64,
    pub core_bottom_temp_degc: f64,
    pub core_top_temp_degc: f64,
    pub core_inlet_temp_degc: f64,
    pub core_outlet_temp_degc: f64,
    pub left_downcomer_upper_temp_degc: f64,
    pub left_downcomer_mid_temp_degc: f64,
    pub left_downcomer_lower_temp_degc: f64,
    pub right_downcomer_upper_temp_degc: f64,
    pub right_downcomer_mid_temp_degc: f64,
    pub right_downcomer_lower_temp_degc: f64,

    // for diagnostics
    /// this displays reactor thermal power in megawatts,
    /// including decay heat
    pub reactor_power_megawatts: f64,
    /// this is decay heat in megawatts 
    pub reactor_decay_heat_megawatts: f64,
    /// this displays reactor keff
    pub keff: f64,
    /// this displays reactivity in dollars 
    pub reactivity_dollars: f64,
    /// this displays xenon feedback in dollars 
    pub xenon135_feedback_dollars: f64,

    // this is important for coupling between prke loop and thermal 
    // hydraulics loop
    pub prke_loop_accumulated_timestep_seconds: f64,
    pub prke_loop_accumulated_heat_removal_kilojoules: f64,

    /// pump pressure settings 
    pub fhr_pri_loop_pump_pressure_kilopascals: f64,
    pub fhr_intermediate_loop_pump_pressure_kilopascals: f64,


    // this is important for timestep monitoring 
    // time diagnostics
    pub prke_simulation_time_seconds: f64,
    pub prke_elapsed_time_seconds: f64,
    pub prke_calc_time_microseconds: f64,
    pub prke_timestep_microseconds: f64,

    pub thermal_hydraulics_simulation_time_seconds: f64,
    pub thermal_hydraulics_calc_time_microseconds: f64,
    pub thermal_hydraulics_timestep_microseconds: f64,


    // mass flowrate 
    // diagnostics for thermal hydraulics loop 
    pub reactor_branch_flowrate_kg_per_s: f64,
    pub downcomer1_branch_flowrate_kg_per_s: f64,
    pub downcomer2_branch_flowrate_kg_per_s: f64,
    pub ihx_branch_flowrate_kg_per_s: f64,
    pub intermediate_loop_clockwise_flow_kg_per_s: f64,


    // temperatures for primary loop diagnostics
    pub pipe_4_temperature_vector_degc: Vec<f64>,
    pub pipe_5_temperature_vector_degc: Vec<f64>,
    pub ihx_shell_6_temperature_vector_degc: Vec<f64>,
    pub pipe_7_temperature_vector_degc: Vec<f64>,
    pub pipe_8_temperature_vector_degc: Vec<f64>,
    pub pri_pump_9_temperature_vector_degc: Vec<f64>,
    pub pipe_10_temperature_vector_degc: Vec<f64>,
    pub pipe_11_temperature_vector_degc: Vec<f64>,

    // temperatures for secondary loop
    pub ihx_tube_6_temperature_vector_degc: Vec<f64>,
    pub pipe_12_temperature_vector_degc: Vec<f64>,
    pub pipe_13_temperature_vector_degc: Vec<f64>,
    /// steam generator (sg) shell side labelled 14,
    /// this is the temperature of 
    /// HITEC that fills the shell side of 
    /// the steam generator
    pub sg_shell_14_temperature_vector_degc: Vec<f64>,
    pub pipe_15_temperature_vector_degc: Vec<f64>,
    pub intrmd_pump_16_temperature_vector_degc: Vec<f64>,
    pub pipe_17_temperature_vector_degc: Vec<f64>,

    // settings for steam generator loop
    pub user_specified_secondary_loop_mass_flowrate_kg_per_s: f64,
    pub user_specified_secondary_loop_pump_outlet_pressure_bar: f64,
    pub user_specified_secondary_loop_ua_watt_per_kelvin:f64,
    pub steam_generator_tube_outlet_temperature_degc: f64,

    // void fractions for secondary loop 
    pub steam_quality_after_condenser: f64,
    pub steam_quality_after_pump: f64,
    pub steam_quality_after_steam_generator_tube_side: f64,
    pub steam_quality_after_turbine: f64,

    // performance metrics for steam cycle 
    pub turbine_power_megawatts: f64,
    pub condenser_duty_megawatts: f64,
}

impl Default for FHRState {
    fn default() -> Self {
        let default_temperature_degc = 500.0;
        FHRState { 
            left_cr_insertion_frac: 0.40,
            right_cr_insertion_frac: 0.40,
            pebble_core_temp_degc: default_temperature_degc,
            pebble_bed_coolant_temp_degc: default_temperature_degc,
            core_bottom_temp_degc: default_temperature_degc,
            core_top_temp_degc: default_temperature_degc,
            core_inlet_temp_degc: default_temperature_degc,
            core_outlet_temp_degc: default_temperature_degc,
            left_downcomer_upper_temp_degc: default_temperature_degc,
            left_downcomer_mid_temp_degc: default_temperature_degc,
            left_downcomer_lower_temp_degc: default_temperature_degc,
            right_downcomer_upper_temp_degc: default_temperature_degc,
            right_downcomer_mid_temp_degc: default_temperature_degc,
            right_downcomer_lower_temp_degc: default_temperature_degc,
            prke_loop_accumulated_timestep_seconds: 0.0,
            prke_loop_accumulated_heat_removal_kilojoules: 0.0,
            reactor_power_megawatts: 0.0,
            keff: 0.0,
            reactivity_dollars: 0.0,
            xenon135_feedback_dollars: 0.0,
            prke_simulation_time_seconds: 0.0,
            prke_elapsed_time_seconds: 0.0,
            prke_calc_time_microseconds: 0.0,
            prke_timestep_microseconds: 0.0,
            reactor_decay_heat_megawatts: 0.0,
            fhr_pri_loop_pump_pressure_kilopascals: 100.0,
            fhr_intermediate_loop_pump_pressure_kilopascals: 100.0,
            thermal_hydraulics_simulation_time_seconds: 0.0,
            thermal_hydraulics_calc_time_microseconds: 0.0,
            thermal_hydraulics_timestep_microseconds: 0.0,
            reactor_branch_flowrate_kg_per_s: 0.0,
            downcomer1_branch_flowrate_kg_per_s: 0.0,
            downcomer2_branch_flowrate_kg_per_s: 0.0,
            ihx_branch_flowrate_kg_per_s: 0.0,
            intermediate_loop_clockwise_flow_kg_per_s: 0.0,
            user_specified_secondary_loop_mass_flowrate_kg_per_s: 50.0,
            user_specified_secondary_loop_pump_outlet_pressure_bar: 1.2,
            user_specified_secondary_loop_ua_watt_per_kelvin: 1.5e5,
            steam_generator_tube_outlet_temperature_degc: 300.0,
            pipe_4_temperature_vector_degc: vec![500.0],
            pipe_5_temperature_vector_degc: vec![500.0],
            ihx_shell_6_temperature_vector_degc: vec![500.0, 500.0],
            pipe_7_temperature_vector_degc: vec![500.0],
            pipe_8_temperature_vector_degc: vec![500.0],
            pri_pump_9_temperature_vector_degc: vec![500.0],
            pipe_10_temperature_vector_degc: vec![500.0],
            pipe_11_temperature_vector_degc: vec![500.0],
            ihx_tube_6_temperature_vector_degc: vec![500.0, 500.0],
            pipe_12_temperature_vector_degc: vec![500.0],
            pipe_13_temperature_vector_degc: vec![500.0],
            sg_shell_14_temperature_vector_degc: vec![500.0],
            pipe_15_temperature_vector_degc: vec![500.0],
            intrmd_pump_16_temperature_vector_degc: vec![500.0],
            pipe_17_temperature_vector_degc: vec![500.0],
            steam_quality_after_condenser: 0.0,
            steam_quality_after_pump: 0.0,
            steam_quality_after_steam_generator_tube_side: 0.0,
            steam_quality_after_turbine: 0.0,
            turbine_power_megawatts: 0.0,
            condenser_duty_megawatts: 0.0,
        }
    }
}

impl FHRState {

    pub fn obtain_average_heat_removal_rate_from_pebble_bed_and_reset_counter(
        &mut self) -> Power {
        let heat_removal_rate_kilowatts = 
            self.prke_loop_accumulated_heat_removal_kilojoules/
            self.prke_loop_accumulated_timestep_seconds;

        self.prke_loop_accumulated_timestep_seconds = 0.0;
        self.prke_loop_accumulated_heat_removal_kilojoules = 0.0;
        return Power::new::<kilowatt>(heat_removal_rate_kilowatts);
    }
}


impl FHRSimulatorApp {
    /// Called once before the first frame.
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        // This is also where you can customize the look and feel of egui using
        // `cc.egui_ctx.set_visuals` and `cc.egui_ctx.set_fonts`.

        //// Load previous app state (if any).
        //// Note that you must enable the `persistence` feature for this to work.
        //if let Some(storage) = cc.storage {
        //    return eframe::get_value(storage, eframe::APP_KEY).unwrap_or_default();
        //}

        let new_fhr_app: FHRSimulatorApp = Default::default();

        let fhr_state_prke_ptr: Arc<Mutex<FHRState>> = 
            new_fhr_app.fhr_state.clone();
        let fhr_state_thermal_hydraulics_ptr: Arc<Mutex<FHRState>> = 
            new_fhr_app.fhr_state.clone();
        // these are pointers/references for plotting reactor power 
        // both the instantaneous state 
        // and page plotting
        let fhr_state_plot_ptr: Arc<Mutex<FHRState>> = 
            new_fhr_app.fhr_state.clone();
        let fhr_page_plot_ptr: Arc<Mutex<PagePlotData>> = 
            new_fhr_app.fhr_simulator_ptr_for_plotting.clone();

        // now spawn a thread to do the kinetics
        //
        thread::spawn(move ||{
            // now I also have a PRKE data which lives inside this loop
            FHRSimulatorApp::calculate_prke_loop(fhr_state_prke_ptr);
        });

        // spawn a thread to do the thermal hydraulics
        thread::spawn(move ||{
            FHRSimulatorApp::calculate_thermal_hydraulics_loop(
                fhr_state_thermal_hydraulics_ptr
            );
            
        });
        // spawn a thread to do the updating of graph plots
        thread::spawn(move ||{
            FHRSimulatorApp::update_plot_from_fhr_state(
                fhr_state_plot_ptr,
                fhr_page_plot_ptr
            );
            
        });

        new_fhr_app
    }


    
}
impl Default for FHRSimulatorApp {
    fn default() -> Self {

        let fhr_state = FHRState::default();
        let fhr_state_ptr = Arc::new(Mutex::new(fhr_state));
        let fhr_plot: PagePlotData = PagePlotData::default();
        let fhr_plot_ptr = Arc::new(Mutex::new(fhr_plot));
        let default_open_panel = Panel::MainPage;

        Self {
            fhr_state: fhr_state_ptr,
            open_panel: default_open_panel,
            fhr_simulator_ptr_for_plotting: fhr_plot_ptr,

        }
    }
}

pub mod app;
