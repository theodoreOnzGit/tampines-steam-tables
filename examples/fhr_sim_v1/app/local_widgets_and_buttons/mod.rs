use crate::FHRSimulatorApp;

use egui::{Pos2, Rect, Ui, Widget};

pub mod reactor_art;

pub mod fhr_reactor_widget;

pub mod pipes;

pub mod pumps;


impl FHRSimulatorApp {
    // places a widget at some area
    pub fn put_widget_with_size_and_centre(
        &mut self, ui: &mut Ui, widget: impl Widget,
        centre_x_pixels: f32,
        centre_y_pixels: f32,
        x_width_pixels: f32,
        y_width_pixels: f32){

        let top_left_x: f32 = centre_x_pixels - 0.5 * x_width_pixels;
        let top_left_y: f32 = centre_y_pixels - 0.5 * y_width_pixels;
        let bottom_right_x: f32 = centre_x_pixels + 0.5 * x_width_pixels;
        let bottom_right_y: f32 = centre_y_pixels + 0.5 * y_width_pixels;

        let rect: Rect = Rect {
            // top left
            min: Pos2 { x: top_left_x, y: top_left_y },
            // bottom right
            max: Pos2 { x: bottom_right_x, y: bottom_right_y },
        };

        ui.put(rect, widget);

    }

    pub fn place_vertical_widget_with_length(
        &mut self, ui: &mut Ui, widget: impl Widget,
        centre_x_pixels: f32,
        centre_y_pixels: f32,
        button_length: f32,
        aspect_ratio: f32,
        ){

        // aspect ratio is length by breadth (longer side by shorter side)
        
        let y_width_pixels = button_length;
        let mut x_width_pixels = button_length/aspect_ratio;

        // min width is 30 px 
        if x_width_pixels < 30.0 {
            x_width_pixels = 30.0;
        }

        self.put_widget_with_size_and_centre(
            ui, 
            widget, 
            centre_x_pixels, 
            centre_y_pixels, 
            x_width_pixels, 
            y_width_pixels);
    }

    pub fn place_horizontal_widget_with_length(
        &mut self, ui: &mut Ui, widget: impl Widget,
        centre_x_pixels: f32,
        centre_y_pixels: f32,
        button_length: f32,
        aspect_ratio: f32,
        ){

        // aspect ratio is length by breadth (longer side by shorter side)
        
        let x_width_pixels = button_length;
        let mut y_width_pixels = button_length/aspect_ratio;
        // min width is 30 px 
        if y_width_pixels < 30.0 {
            y_width_pixels = 30.0;
        }

        self.put_widget_with_size_and_centre(
            ui, 
            widget, 
            centre_x_pixels, 
            centre_y_pixels, 
            x_width_pixels, 
            y_width_pixels);
    }

    
}



pub fn new_temp_sensitive_button_blue_red(
    min_temp_degc: f32, 
    max_temp_degc: f32,
    button_temp_degc: f32,
    name: &str,
) -> egui::Button {

    let hotness: f32 = 
        (button_temp_degc - min_temp_degc)/(max_temp_degc- min_temp_degc);

    let colour_temp = hot_to_cold_colour_mark_1(hotness);
    let temp_sensitive_button = egui::Button::new(name)
        .fill(colour_temp);

    temp_sensitive_button

}

pub fn new_temp_sensitive_button_black_red(
    min_temp_degc: f32, 
    max_temp_degc: f32,
    button_temp_degc: f32,
    name: &str,
) -> egui::Button {

    let hotness: f32 = 
        (button_temp_degc - min_temp_degc)/(max_temp_degc- min_temp_degc);

    let colour_temp = hot_to_cold_colour_mark_2(hotness);
    let temp_sensitive_button = egui::Button::new(name)
        .fill(colour_temp);

    temp_sensitive_button

}




/// contains traits to help objects in the simulation 
/// visually connect to one another
pub mod simulator_trait;

/// hot to cold colours 
pub mod hot_to_cold_colour_functions;
pub use hot_to_cold_colour_functions::*;

/// steam quality colours  (blue for 0, white for 1)
pub mod steam_quality_colour_functions;

/// turbine art 
pub mod turbine_widget;
