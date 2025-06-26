use std::{sync::{Arc, Mutex}, thread};

use uom::si::{f64::*, power::kilowatt};

use crate::app::{graph_data::PagePlotData, panel_enum::Panel};


/// this represents the first iteration 
/// of the fhr simulator
///
/// basically one can do a FHR loop 
/// with a permenantly steady state steam cycle
/// the latter uses the tampines-steam-tables
#[cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release
fn main(){

    fhr_simulator_v1().unwrap();


}
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
#[derive(Clone,Copy, Debug)]
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


    // diagnostics for thermal hydraulics loop 
    pub reactor_branch_flowrate_kg_per_s: f64,
    pub downcomer1_branch_flowrate_kg_per_s: f64,
    pub downcomer2_branch_flowrate_kg_per_s: f64,
    pub ihx_branch_flowrate_kg_per_s: f64,
    pub intermediate_loop_clockwise_flow_kg_per_s: f64,
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
