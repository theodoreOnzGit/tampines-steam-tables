use std::ops::Deref;

use egui::Ui;

use crate::{FHRSimulatorApp, FHRState};

use super::local_widgets_and_buttons::new_temp_sensitive_button;

impl FHRSimulatorApp {

    pub(crate) fn side_panel(&mut self, ui: &mut Ui){

        egui::ScrollArea::both()
            .scroll_bar_visibility(egui::scroll_area::ScrollBarVisibility::AlwaysVisible)
            .drag_to_scroll(true)
            .show(ui, |ui| {

                ui.heading("Fluoride Salt Cooled High Temperature Reactor (FHR) Controls");
                let mut fhr_state_ptr = self.fhr_state.lock().unwrap();

                let left_cr_slider = egui::Slider::new(
                    &mut fhr_state_ptr.left_cr_insertion_frac, 
                    0.0000..=1.0)
                    .logarithmic(false)
                    .text("Left Control Rod insertion Fraction")
                    .drag_value_speed(0.001);

                ui.add(left_cr_slider);

                let right_cr_slider = egui::Slider::new(
                    &mut fhr_state_ptr.right_cr_insertion_frac, 
                    0.0000..=1.0)
                    .logarithmic(false)
                    .text("Right Control Rod insertion Fraction")
                    .drag_value_speed(0.001);

                ui.add(right_cr_slider);

                let pri_loop_pump_slider = egui::Slider::new(
                    &mut fhr_state_ptr.fhr_pri_loop_pump_pressure_kilopascals, 
                    0.0..=200.0)
                    .logarithmic(false)
                    .text("Primary Loop Pump Pressure (kPa)")
                    .drag_value_speed(0.001);

                ui.add(pri_loop_pump_slider);
                let intermediate_loop_pump_slider = egui::Slider::new(
                    &mut fhr_state_ptr.fhr_intermediate_loop_pump_pressure_kilopascals, 
                    0.0..=200.0)
                    .logarithmic(false)
                    .text("Intermediate Loop Pump Pressure (kPa)")
                    .drag_value_speed(0.001);

                ui.add(intermediate_loop_pump_slider);


                // cloning the entire fhr state for diagnostics
                let fhr_state_clone: FHRState = fhr_state_ptr.deref().clone();
                //
                drop(fhr_state_ptr);

                ui.separator();
                ui.heading("FHR Diagnostics");

                let pebble_core_temp_degc = 
                    fhr_state_clone.pebble_core_temp_degc;

                let pebble_bed_coolant_temp_degc = 
                    fhr_state_clone.pebble_bed_coolant_temp_degc;

                // need pebble bed power and/or heat removal
                // and keff
                let keff = fhr_state_clone.keff;
                let reactor_power_megawatts = 
                    fhr_state_clone.reactor_power_megawatts;

                ui.label("Reactor Power (MW-thermal):");
                ui.label(((1000.0*reactor_power_megawatts).round() / 1000.0).to_string());

                ui.label("Fuel Temperature Pebble Core/TRISO (deg C):");
                ui.label(((10.0*pebble_core_temp_degc).round() / 10.0).to_string());
                ui.label("Pebble Bed Coolant Temp (deg C):");
                ui.label(((10.0*pebble_bed_coolant_temp_degc).round() / 10.0).to_string());
                ui.label("k_eff");
                ui.label(((1.0e6*keff).round() / 1.0e6).to_string());

                let reactivity_dollars = fhr_state_clone.reactivity_dollars;
                ui.label("Reactivity ($ dollars)");
                ui.label(((1.0e3*reactivity_dollars).round() / 1.0e3).to_string());

                let xe135_feedback_dollars = fhr_state_clone.xenon135_feedback_dollars;
                ui.label("Xe135 feedback ($ dollars)");
                ui.label(((1.0e3*xe135_feedback_dollars).round() / 1.0e3).to_string());

                // then temperature scale 

                ui.separator();
                ui.heading("Temperature Scale");
                ui.heading("Colour to Temperature Legend");

                // now I need colour legend
                let min_temp_degc = 450.0;
                let max_temp_degc = 1000.0;
                // max temp
                let button_temp_degc = max_temp_degc;
                let max_temp_string: String = 
                    button_temp_degc.to_string()+" degC or more";
                let max_temp = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &max_temp_string
                );

                ui.add(max_temp);
                // 950.0
                let button_temp_degc = 950.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_950_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_950_degc);
                // 900.0
                let button_temp_degc = 900.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_900_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_900_degc);
                // 850.0
                let button_temp_degc = 850.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_850_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_850_degc);
                // 800.0
                let button_temp_degc = 800.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_800_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_800_degc);
                // 750.0
                let button_temp_degc = 750.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_750_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_750_degc);

                // 700.0
                let button_temp_degc = 700.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_700_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_700_degc);

                // 650.0
                let button_temp_degc = 650.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_650_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_650_degc);
                // 600.0
                let button_temp_degc = 600.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_600_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_600_degc);
                // 550.0
                let button_temp_degc = 550.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_550_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_550_degc);
                // 500.0
                let button_temp_degc = 500.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_500_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_500_degc);
                // 450.0
                let button_temp_degc = 450.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_450_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_450_degc);

                // 400.0
                let button_temp_degc = 400.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_400_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_400_degc);

                // 350.0
                let button_temp_degc = 350.0;
                let button_temp_string: String = 
                    button_temp_degc.to_string()+" degrees celsius";
                let temp_350_degc = new_temp_sensitive_button(
                    min_temp_degc, 
                    max_temp_degc, 
                    button_temp_degc, 
                    &button_temp_string
                );
                ui.add(temp_350_degc);

                // time diagnostics 
                ui.separator();
                ui.heading("Timestep Diagnostics");
                ui.separator();

                let prke_elapsed_time_seconds = fhr_state_clone.prke_elapsed_time_seconds;

                ui.label("Elapsed Time Seconds");
                ui.label(((1000.0*prke_elapsed_time_seconds).round() / 1000.0).to_string());


                let prke_simulation_time_seconds = fhr_state_clone.prke_simulation_time_seconds;

                ui.label("PRKE Simulation Time Seconds");
                ui.label(((1000.0*prke_simulation_time_seconds).round() / 1000.0).to_string());


                let prke_timestep_microseconds = fhr_state_clone.prke_timestep_microseconds;

                ui.label("PRKE Timestep Microseconds");
                ui.label(((1000.0*prke_timestep_microseconds).round() / 1000.0).to_string());

                let prke_calc_time_microseconds = fhr_state_clone.prke_calc_time_microseconds;

                ui.label("PRKE Calculation time per timestep Microseconds");
                ui.label(((1000.0*prke_calc_time_microseconds).round() / 1000.0).to_string());


                ui.separator();
                ui.separator();


                let thermal_hydraulics_simulation_time_seconds = fhr_state_clone.thermal_hydraulics_simulation_time_seconds;

                ui.label("Thermal Hydraulics Simulation Time Seconds");
                ui.label(((1000.0*thermal_hydraulics_simulation_time_seconds).round() / 1000.0).to_string());

                let thermal_hydraulics_timestep_microseconds = fhr_state_clone.thermal_hydraulics_timestep_microseconds;

                ui.label("Thermal Hydraulics Timestep Microseconds");
                ui.label(((1000.0*thermal_hydraulics_timestep_microseconds).round() / 1000.0).to_string());

                let thermal_hydraulics_calc_time_microseconds = fhr_state_clone.thermal_hydraulics_calc_time_microseconds;

                ui.label("Thermal Hydraulics Calculation time per timestep Microseconds");
                ui.label(((1000.0*thermal_hydraulics_calc_time_microseconds).round() / 1000.0).to_string());



                ui.separator();

                // flowrate diagnostics
                ui.separator();
                ui.heading("Mass Flowrate Diagnostics");
                ui.separator();


                let ihx_br_flowrate_kg_per_s = -fhr_state_clone.ihx_branch_flowrate_kg_per_s;

                ui.label("Intermediate Heat Exchanger Branch Flowrate downwards through IHX (kg/s)");
                ui.label(((1000.0*ihx_br_flowrate_kg_per_s).round() / 1000.0).to_string());

                let reactor_branch_flowrate_kg_per_s = fhr_state_clone.reactor_branch_flowrate_kg_per_s;

                ui.label("Reactor Branch Flowrate (upwards) (kg/s)");
                ui.label(((1000.0*reactor_branch_flowrate_kg_per_s).round() / 1000.0).to_string());

                let downcomer1_branch_flowrate_kg_per_s = fhr_state_clone.downcomer1_branch_flowrate_kg_per_s;

                ui.label("Left Downcomer Branch Flowrate (upwards) (kg/s)");
                ui.label(((1000.0*downcomer1_branch_flowrate_kg_per_s).round() / 1000.0).to_string());


                let downcomer2_branch_flowrate_kg_per_s = fhr_state_clone.downcomer2_branch_flowrate_kg_per_s;

                ui.label("Right Downcomer Branch Flowrate (upwards) (kg/s)");
                ui.label(((1000.0*downcomer2_branch_flowrate_kg_per_s).round() / 1000.0).to_string());

                let intermediate_loop_clockwise_flowrate_kg_per_s = fhr_state_clone.intermediate_loop_clockwise_flow_kg_per_s;

                ui.label("Intermediate Loop Flowrate (clockwise) (kg/s)");
                ui.label(((1000.0*intermediate_loop_clockwise_flowrate_kg_per_s).round() / 1000.0).to_string());
                // then acknowledgements/citing

            });


    }
}
