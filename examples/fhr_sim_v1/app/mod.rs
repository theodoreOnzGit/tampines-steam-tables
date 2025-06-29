use std::time::Duration;

use egui::{vec2, Pos2, Rect, Vec2};
use local_widgets_and_buttons::{fhr_reactor_widget::FHRReactorWidget, pipes::SinglePipe};
use uom::si::f64::*;
use uom::si::thermodynamic_temperature::degree_celsius;

use crate::{FHRSimulatorApp, FHRState};
use crate::Panel;

pub mod prke_backend;
pub mod thermal_hydraulics_backend;

impl eframe::App for FHRSimulatorApp {
    /// Called by the frame work to save state before shutdown.
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, eframe::APP_KEY, self);
    }

    /// Called each time the UI needs repainting, which may be many times per second.
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Put your widgets into a `SidePanel`, `TopPanel`, `CentralPanel`, `Window` or `Area`.
        // For inspiration and more examples, go to https://emilk.github.io/egui



        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.heading("FHR Educational Simulator v1");
            ui.separator();
            egui::menu::bar(ui, |ui| {

                egui::widgets::global_theme_preference_buttons(ui);
            });
            // allow user to select which panel is open
            ui.horizontal( 
                |ui| {
                    ui.selectable_value(&mut self.open_panel, Panel::MainPage, "Main Page"); 
                    ui.selectable_value(&mut self.open_panel, Panel::ReactorPowerGraphs, "Power Diagnostics"); 
                    ui.selectable_value(&mut self.open_panel, Panel::PoisonGraphs, "Reactor Poison Diagnostics"); 
            }
            );
            ui.separator();
        });

        egui::SidePanel::right("Supplementary Info").show(ctx, |ui|{
            self.side_panel(ui);



        });

        egui::CentralPanel::default().show(ctx, |ui| {

            ui.separator();
            egui::ScrollArea::both()
                .scroll_bar_visibility(egui::scroll_area::ScrollBarVisibility::AlwaysVisible)
                .drag_to_scroll(true)
                .show(ui, |ui| {


                    match self.open_panel {
                        Panel::MainPage => self.main_page(ui),
                        Panel::ReactorPowerGraphs => self.reactor_power_page_graph(ui),
                        Panel::PoisonGraphs => {},
                    }
                });













        });

        egui::TopBottomPanel::bottom("github").show(ctx, |ui|{

            ui.with_layout(egui::Layout::bottom_up(egui::Align::LEFT), |ui| {
                powered_by_egui_and_eframe(ui);
                egui::warn_if_debug_build(ui);
            });

        });

        


        ctx.request_repaint_after(Duration::from_millis(50));

        // adding the return here because there are too many closing 
        // parantheses
        // just demarcates the end
        return ();
    }

}

impl FHRSimulatorApp {
    pub fn main_page(&mut self, ui: &mut egui::Ui,){

        // for painting widgets
        // https://github.com/emilk/egui/blob/master/crates/egui_demo_lib/src/demo/misc_demo_window.rs
        //
        // the main thing is the painter class:
        // https://docs.rs/egui/latest/egui/struct.Painter.html
        //
        // here you can paint circles and rectangles 
        // images, line segments etc.
        // obtain lock first 

        // quickly clone the fhr state and drop ptr asap 
        // just to read
        let fhr_state_ptr = self.fhr_state.lock().unwrap();
        let fhr_state_clone: FHRState = fhr_state_ptr.clone();
        drop(fhr_state_ptr);

        let left_control_rod_insertion_frac 
            = fhr_state_clone.left_cr_insertion_frac;
        let right_control_rod_insertion_frac 
            = fhr_state_clone.right_cr_insertion_frac;



        let ui_rectangle: Rect = ui.min_rect();

        // this gives coordinates of top and left of the ui
        // for relative placement
        let left_most_side = ui_rectangle.left();
        let top_most_side = ui_rectangle.top();

        let reactor_offset_x: f32 = 100.0;
        let reactor_offset_y: f32 = 200.0;
        let reactor_x_width_px: f32 = 150.0 * 1.5;
        let reactor_y_height_px: f32 = 700.0 * 1.5;


        let reactor_rect_top_left: Pos2 = 
            Pos2 { 
                x: left_most_side + reactor_offset_x, 
                y: top_most_side + reactor_offset_y
            };
        let reactor_rect_bottom_right: Pos2 = 
            Pos2 { 
                x: reactor_rect_top_left.x + reactor_x_width_px, 
                y: reactor_rect_top_left.y + reactor_y_height_px
            };
        let reactor_rectangle: egui::Rect =
            egui::Rect{
                min: reactor_rect_top_left,
                max: reactor_rect_bottom_right,
            };


        let fhr_size = 
            vec2(reactor_rectangle.width(), reactor_rectangle.height());

        let min_temp = ThermodynamicTemperature::new::<degree_celsius>(450.0);
        let max_temp = ThermodynamicTemperature::new::<degree_celsius>(750.0);
        let pebble_core_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.pebble_core_temp_degc
        );
        let pebble_bed_coolant_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.pebble_bed_coolant_temp_degc
        );
        let core_bottom_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.core_bottom_temp_degc
        );
        let core_top_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.core_top_temp_degc
        );
        let core_inlet_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.core_inlet_temp_degc
        );
        let core_outlet_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.core_outlet_temp_degc
        );
        let left_downcomer_upper_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.left_downcomer_upper_temp_degc
        );
        let left_downcomer_mid_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.left_downcomer_mid_temp_degc
        );
        let left_downcomer_lower_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.left_downcomer_lower_temp_degc
        );
        let right_downcomer_upper_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.right_downcomer_upper_temp_degc
        );
        let right_downcomer_mid_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.right_downcomer_mid_temp_degc
        );
        let right_downcomer_lower_temp = ThermodynamicTemperature::new::<degree_celsius>(
            fhr_state_clone.right_downcomer_lower_temp_degc
        );

        // pri loop
        let pipe_4_temperature_vector_degc = 
            fhr_state_clone.pipe_4_temperature_vector_degc;
        let pipe_5_temperature_vector_degc = 
            fhr_state_clone.pipe_5_temperature_vector_degc;

        // ihx sthe, "dangerous" temperatures are here 
        let ihx_shell_6_temperature_vector_degc = 
            fhr_state_clone.ihx_shell_6_temperature_vector_degc;
        let ihx_tube_6_temperature_vector_degc = 
            fhr_state_clone.ihx_tube_6_temperature_vector_degc;

        let pipe_7_temperature_vector_degc = 
            fhr_state_clone.pipe_7_temperature_vector_degc;
        let pipe_8_temperature_vector_degc = 
            fhr_state_clone.pipe_8_temperature_vector_degc;
        let pump_9_temperature_vector_degc = 
            fhr_state_clone.pri_pump_9_temperature_vector_degc;
        let pipe_10_temperature_vector_degc = 
            fhr_state_clone.pipe_10_temperature_vector_degc;
        let pipe_11_temperature_vector_degc = 
            fhr_state_clone.pipe_11_temperature_vector_degc;

        // intermediate loop
        let pipe_12_temperature_vector_degc = 
            fhr_state_clone.pipe_12_temperature_vector_degc;
        let pipe_13_temperature_vector_degc = 
            fhr_state_clone.pipe_13_temperature_vector_degc;
        let sg_shell_14_temperature_vector_degc = 
            fhr_state_clone.sg_shell_14_temperature_vector_degc;
        let pipe_15_temperature_vector_degc = 
            fhr_state_clone.pipe_15_temperature_vector_degc;
        let intrmd_pump_16_temperature_vector_degc = 
            fhr_state_clone.intrmd_pump_16_temperature_vector_degc;
        let pipe_17_temperature_vector_degc = 
            fhr_state_clone.pipe_17_temperature_vector_degc;


        let mut fhr_widget = FHRReactorWidget::new(
            fhr_size,
            min_temp,
            max_temp,
            pebble_core_temp,
            pebble_bed_coolant_temp,
            core_bottom_temp,
            core_top_temp,
            core_inlet_temp,
            core_outlet_temp,
            left_downcomer_upper_temp,
            left_downcomer_mid_temp,
            left_downcomer_lower_temp,
            right_downcomer_upper_temp,
            right_downcomer_mid_temp,
            right_downcomer_lower_temp,
        );
        fhr_widget.set_left_cr_frac(left_control_rod_insertion_frac);
        fhr_widget.set_right_cr_frac(right_control_rod_insertion_frac);

        ui.put(reactor_rectangle, fhr_widget);
        
        fn average_temp(temp_vec_degc: &Vec<f64>) -> ThermodynamicTemperature {

            let pipe_temp_sum: f64 = temp_vec_degc.iter().sum();
            let pipe_temp: ThermodynamicTemperature = 
                ThermodynamicTemperature::new::<degree_celsius>(
                    pipe_temp_sum/temp_vec_degc.len() as f64
                );

            return pipe_temp;
        }
        // these are reactor lengthscales to pay attention to 
        // 
        let reactor_height: f32 = reactor_rectangle.height() * 0.56;
        let reactor_width: f32 = reactor_rectangle.width();

        let pipe_11_temp: ThermodynamicTemperature = 
            average_temp(&pipe_11_temperature_vector_degc);
        // this represents the pipe coordinate change from the start 
        // point
        let pipe_11_coordinate_chg_as_percentage_of_reactor = 
            vec2(0.0, 30.0);

        // the start point of the pipe 11 
        // is: x coordinate (boundy by left and right of reactor rectangle)
        // y coordinate: reactor_rectangle bottom minus some fixed height 
        // based on scale
        let pipe_11_start = 
            vec2(
                0.5 * reactor_rectangle.left() + 0.5 * reactor_rectangle.right(),
                reactor_rectangle.bottom() - reactor_rectangle.height() * 0.28,
            );

        let pipe_10_start = 
            vec2(
                pipe_11_start.x + pipe_11_coordinate_chg_as_percentage_of_reactor.x/100.0 * reactor_width,
                pipe_11_start.y + pipe_11_coordinate_chg_as_percentage_of_reactor.y/100.0 * reactor_height,
            );


        let pipe_11_rect = 
            egui::Rect {
                min: Pos2 { x: 0.0, y: 0.0 } + pipe_11_start,
                max: Pos2 { x: 0.0, y: 0.0 } + pipe_10_start,
            };


        let pipe_11_coordinate_chg = 
            vec2(
                pipe_11_coordinate_chg_as_percentage_of_reactor.x/100.0 * reactor_width,
                pipe_11_coordinate_chg_as_percentage_of_reactor.y/100.0 * reactor_height,
            );

        let pipe_11_widget = SinglePipe::new(
            pipe_11_coordinate_chg, 
            min_temp, 
            max_temp, 
            pipe_11_temp
        );
        ui.put(pipe_11_rect, pipe_11_widget);


        fn create_pipe_widget (
            pipe_temp_vec_degc: &Vec<f64>,
            start_point: Vec2, 
            pipe_position_change_as_percentage_of_reactor: Vec2,
            ui: &mut egui::Ui,
            reactor_width: f32,
            reactor_height: f32) -> Vec2 {


                // the start point of the pipe 11 
                // is: x coordinate (boundy by left and right of reactor rectangle)
                // y coordinate: reactor_rectangle bottom minus some fixed height 
                // based on scale

                let pipe_temp = 
                    average_temp(pipe_temp_vec_degc);
                let end_point = 
                    vec2(
                        start_point.x + pipe_position_change_as_percentage_of_reactor.x/100.0 * reactor_width,
                        start_point.y + pipe_position_change_as_percentage_of_reactor.y/100.0 * reactor_height,
                    );


                // now in the case that end point is 
                // higher in x and y position than the start point:
                let mut pipe_rect = 
                    egui::Rect {
                        min: Pos2 { x: 0.0, y: 0.0 } + start_point,
                        max: Pos2 { x: 0.0, y: 0.0 } + end_point,
                    };
                // this may not always be the case, and the min 
                // must be the top left corner, the end point must 
                // be the bottom right corner
                // 
                // so end point is bottom right if 
                // its x and y coordinate are more than the start 
                // point 
                let end_point_is_bottom_right: bool = 
                    end_point.x > start_point.x && 
                    end_point.y > start_point.y;

                if end_point_is_bottom_right {

                    // do nothing in this case, no need to modify the 
                    // pipe rectangle
                } else {
                    // we need to find the bottom right first 

                    let mut top_left_x_coord: f32 = start_point.x;
                    let mut bottom_right_x_coord: f32 = end_point.x;
                    // 
                    // if the start point happens to be at the bottom right 
                    // so to speak,
                    if end_point.x < start_point.x {
                        top_left_x_coord = end_point.x;
                        bottom_right_x_coord = start_point.x;
                    };

                    let mut top_left_y_coord: f32 = start_point.y;
                    let mut bottom_right_y_coord: f32 = end_point.y;
                    // 
                    // if the start point happens to be at the bottom right 
                    // so to speak,
                    if end_point.y < start_point.y {
                        top_left_y_coord = end_point.y;
                        bottom_right_y_coord = start_point.y;
                    };

                    let top_left = 
                        vec2(
                            top_left_x_coord,
                            top_left_y_coord,
                        );
                    let bottom_right = 
                        vec2(
                            bottom_right_x_coord,
                            bottom_right_y_coord,
                        );

                    pipe_rect = 
                        egui::Rect {
                            min: Pos2 { x: 0.0, y: 0.0 } + top_left,
                            max: Pos2 { x: 0.0, y: 0.0 } + bottom_right,
                        };


                }

                let pipe_coordinate_chg = 
                    vec2(
                        pipe_position_change_as_percentage_of_reactor.x/100.0 * reactor_width,
                        pipe_position_change_as_percentage_of_reactor.y/100.0 * reactor_height,
                    );

                let min_temp = ThermodynamicTemperature::new::<degree_celsius>(450.0);
                let max_temp = ThermodynamicTemperature::new::<degree_celsius>(750.0);

                let pipe_widget = SinglePipe::new(
                    pipe_coordinate_chg, 
                    min_temp, 
                    max_temp, 
                    pipe_temp
                );
                ui.put(pipe_rect, pipe_widget);

                return end_point;
            }



        // now let's create pipe_10 
        let pipe_10_coordinate_chg_percentage = 
            vec2(100.0, 0.0);
        let pump_9_start_point = create_pipe_widget(
            &pipe_10_temperature_vector_degc,
            pipe_10_start,
            pipe_10_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        // pump 9 
        let pump_9_coordinate_chg_percentage = 
            vec2(100.0, 0.0);
        let pipe_8_start_point = create_pipe_widget(
            &pump_9_temperature_vector_degc,
            pump_9_start_point,
            pump_9_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );
        ui.separator();
        // pipe 8
        let pipe_8_coordinate_chg_percentage = 
            vec2(0.0, -30.0);
        let pipe_7_start_point = create_pipe_widget(
            &pipe_8_temperature_vector_degc,
            pipe_8_start_point,
            pipe_8_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );
        // pipe 7
        let pipe_7_coordinate_chg_percentage = 
            vec2(0.0, -100.0);
        let ihx_shell_6_start_point = create_pipe_widget(
            &pipe_7_temperature_vector_degc,
            pipe_7_start_point,
            pipe_7_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );
        // ihx shell
        let ihx_shell_6_coordinate_chg_percentage = 
            vec2(0.0, -30.0);
        let pipe_5_start_point = create_pipe_widget(
            &ihx_shell_6_temperature_vector_degc,
            ihx_shell_6_start_point,
            ihx_shell_6_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );
        // pipe_5
        let pipe_5_coordinate_chg_percentage = 
            vec2(-200.0, -0.0);
        let pipe_4_start_point = create_pipe_widget(
            &pipe_5_temperature_vector_degc,
            pipe_5_start_point,
            pipe_5_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );
        // pipe_4
        let pipe_4_coordinate_chg_percentage = 
            vec2(-0.0, 60.0);
        let _pipe_3_start_point = create_pipe_widget(
            &pipe_4_temperature_vector_degc,
            pipe_4_start_point,
            pipe_4_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        // intermediate loop

        let ihx_tube_6_start_point = 
            ihx_shell_6_start_point +
            vec2(reactor_width * 0.1, 0.0);

        let ihx_tube_6_coordinate_chg_percentage = 
            vec2(0.0, -30.0);
        let ihx_tube_6b_start_point = create_pipe_widget(
            &ihx_tube_6_temperature_vector_degc,
            ihx_tube_6_start_point,
            ihx_tube_6_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );


        // this is tube to make it curl back from heat exchanger
        let ihx_tube_6_coordinate_chg_percentage = 
            vec2(30.0, 0.0);
        let pipe_17_start_point = create_pipe_widget(
            &ihx_tube_6_temperature_vector_degc,
            ihx_tube_6b_start_point,
            ihx_tube_6_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        let ihx_tube_6a_end_point = 
            ihx_tube_6_start_point;

        let ihx_tube_6a_start_point = create_pipe_widget(
            &ihx_tube_6_temperature_vector_degc,
            ihx_tube_6a_end_point,
            ihx_tube_6_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        // pipe 17
        let pipe_17_coordinate_chg_percentage = 
            vec2(150.0, 0.0);
        let _pipe_17_end_point = create_pipe_widget(
            &pipe_17_temperature_vector_degc,
            pipe_17_start_point,
            pipe_17_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        // pipe 12 
        let pipe_12_end_point = ihx_tube_6a_start_point;

        let pipe_12_coordinate_chg_percentage = 
            vec2(0.0, -130.0);
        let pipe_12_start_point = pipe_12_end_point 
            - vec2(
                reactor_width * pipe_12_coordinate_chg_percentage.x/100.0, 
                reactor_height * pipe_12_coordinate_chg_percentage.y/100.0,
            );

        let _pipe_12_end_point = create_pipe_widget(
            &pipe_12_temperature_vector_degc,
            pipe_12_start_point,
            pipe_12_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        // steam generator branch
        //
        // pump 16
        let pump_16_start_point = pipe_12_start_point;

        let pump_16_coordinate_chg_percentage = 
            vec2(75.0, 0.0);
        let pipe_15_start_point = create_pipe_widget(
            &intrmd_pump_16_temperature_vector_degc,
            pump_16_start_point,
            pump_16_coordinate_chg_percentage,
            ui,
            reactor_width,
            reactor_height,
        );

        // pipe 15
        let pipe_15_coordinate_chg_percentage = 
            vec2(75.0, 0.0);
        let sg_shell_14_start_point = create_pipe_widget(
            &pipe_15_temperature_vector_degc, 
            pipe_15_start_point, 
            pipe_15_coordinate_chg_percentage, 
            ui, reactor_width, reactor_height);

        // steam generator (sg) shell side 14
        let sg_shell_14_coordinate_chg_percentage = 
            vec2(0.0, -30.0);
        let pipe_13_start_point = create_pipe_widget(
            &sg_shell_14_temperature_vector_degc, 
            sg_shell_14_start_point, 
            sg_shell_14_coordinate_chg_percentage, 
            ui, reactor_width, reactor_height);

        // pipe 13 
        let pipe_13_coordinate_chg_percentage = 
            vec2(0.0, -130.0);
        let pipe_13_start_point = create_pipe_widget(
            &pipe_13_temperature_vector_degc, 
            pipe_13_start_point, 
            pipe_13_coordinate_chg_percentage, 
            ui, reactor_width, reactor_height);


    }

}

fn powered_by_egui_and_eframe(ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        ui.label("Powered by ");
        ui.hyperlink_to("egui", "https://github.com/emilk/egui");
        ui.label(" and ");
        ui.hyperlink_to(
            "eframe",
            "https://github.com/emilk/egui/tree/master/crates/eframe",
        );
        ui.label(".");
    });
}



pub mod local_widgets_and_buttons;

pub mod side_panel;

pub mod panel_enum;

pub mod graph_data;

pub mod graph_pages;
