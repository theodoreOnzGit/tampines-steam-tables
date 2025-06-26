use egui::Ui;
use egui_plot::{Legend, Line, Plot, PlotPoints};

use crate::{app::graph_data::PagePlotData, FHRSimulatorApp};


impl FHRSimulatorApp {
    pub fn reactor_power_page_graph(&mut self, ui: &mut Ui){

        ui.horizontal(|ui| {
            ui.label("Reactor Power Page");
            if ui.button("Update CSV Data").clicked(){
                //// spawn a new window with csv data
                //let latest_ciet_plot_data: PagePlotData = 
                //    self.ciet_plot_data_mutex_ptr_for_parallel_data_transfer
                //    .lock().unwrap().clone();

                //self.ciet_plot_data = latest_ciet_plot_data;

            };
        });
        ui.separator();
        ui.separator();
        egui::ScrollArea::both().show(ui, |ui| {
            let mut reactor_power_plot = Plot::new("heater inlet and outlet temp degC").legend(Legend::default());

            // sets the aspect for plot 
            reactor_power_plot = reactor_power_plot.width(800.0);
            reactor_power_plot = reactor_power_plot.view_aspect(16.0/9.0);

            reactor_power_plot = reactor_power_plot.x_axis_label(
                "time (seconds), current time (seconds): ".to_owned() 
            );
            reactor_power_plot = reactor_power_plot.y_axis_label(
                "Reactor Power (Megawatts)".to_owned());
            let latest_fhr_plot_data: PagePlotData = 
                self.fhr_simulator_ptr_for_plotting.lock().unwrap().clone();


            // let's make the time and reactor power vector
            // that is with and without decay heat
            let time_reactor_power_vec: Vec<[f64;2]> = 
                latest_fhr_plot_data.get_reactor_power_megawatts_vs_time_secs_vec();

            let time_reactor_power_no_decay_heat_vec: Vec<[f64;2]> = 
                latest_fhr_plot_data.get_reactor_power_no_decay_heat_megawatts_vs_time_secs_vec();

            ui.heading("Reactor Power vs Time");
            reactor_power_plot.show(ui, |plot_ui| {
                plot_ui.line(Line::new(PlotPoints::from(
                            time_reactor_power_vec.clone()
                )).name("Reactor power with decay heat Megawatts"));
                plot_ui.line(Line::new(PlotPoints::from(
                            time_reactor_power_no_decay_heat_vec.clone()
                )).name("Reactor power without decay heat Megawatts"));
                //plot_ui.line(Line::new(PlotPoints::from(
                //            time_simulated_reactor_feedback_outlet_temp_vec.clone()
                //)).name("simulated reactivity bt12 (heater outlet) temperature deg C"));
            });



            //self.citation_disclaimer_and_acknowledgements(ui);

        });

    }
}
