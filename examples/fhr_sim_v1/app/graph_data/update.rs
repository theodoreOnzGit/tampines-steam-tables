use std::{ops::DerefMut, sync::{Arc, Mutex}, time::Duration};


use uom::si::{f64::*, power::megawatt, time::second};

use crate::{FHRSimulatorApp, FHRState};

use super::PagePlotData;

impl FHRSimulatorApp {
    pub fn update_plot_from_fhr_state(
        fhr_state_ptr: Arc<Mutex<FHRState>>,
        fhr_simulator_ptr_for_plotting: Arc<Mutex<PagePlotData>>,
    ){
        loop {
            // first lock fhr state ptr and clone 
            let local_fhr_state: FHRState = 
                fhr_state_ptr.lock().unwrap().clone();

            // get the current plot object
            let mut local_fhr_plot: PagePlotData = 
                fhr_simulator_ptr_for_plotting.lock().unwrap().clone();

            let simulation_time: Time = Time::new::<second>(
                local_fhr_state.prke_simulation_time_seconds);
            {

                let reactor_power_with_decay_heat: Power = 
                    Power::new::<megawatt>(
                        local_fhr_state.reactor_power_megawatts);
                let reactor_decay_heat: Power = 
                    Power::new::<megawatt>(
                        local_fhr_state.reactor_decay_heat_megawatts);

                let reactor_power_without_decay_heat = 
                    reactor_power_with_decay_heat - 
                    reactor_decay_heat.abs();

                local_fhr_plot.insert_reactor_power_data(
                    simulation_time, 
                    reactor_power_with_decay_heat, 
                    reactor_power_without_decay_heat);



            }


            // update after everything
            let data_record_interval_seconds = 
                local_fhr_plot.graph_data_record_interval_seconds;

            let data_record_interval_ms = 
                data_record_interval_seconds * 1000.0 ;

            // update the plot
            *fhr_simulator_ptr_for_plotting.lock().unwrap().deref_mut()
                = local_fhr_plot;

            // historian records every 100ms 
            // or rather, how often the user decides (in intervals)
            std::thread::sleep(Duration::from_millis(data_record_interval_ms as u64));
        }

    }
}
