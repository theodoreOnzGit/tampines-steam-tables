use egui::{vec2, Sense, Stroke, Vec2, Widget};
use uom::si::{f64::*, thermodynamic_temperature::degree_celsius};


use crate::app::local_widgets_and_buttons::{hot_to_cold_colour_mark_2, steam_quality_colour_functions::steam_quality_colour_mark_1};

use super::hot_to_cold_colour_mark_1;

pub struct SinglePipeColourBlueRedTempSensitive {
    size: Vec2,
    vector: Vec2,
    min_temp: ThermodynamicTemperature,
    max_temp: ThermodynamicTemperature,
    temp: ThermodynamicTemperature,
}

impl SinglePipeColourBlueRedTempSensitive {

    /// still need to correct for minimum size
    ///
    /// note that you will put in a vector from start to end 
    /// like how many pixels in x,y direction you want to go 
    /// then the pipe will autosize everything
    pub fn new(vector: Vec2,
        min_temp: ThermodynamicTemperature,
        max_temp: ThermodynamicTemperature,
        temp: ThermodynamicTemperature,) -> Self {

        let min_width = 20.0;

        // now the size here
        let mut size = vector;

        if size.x <= min_width {
            size.x = min_width
        };

        if size.y <= min_width {
            size.y = min_width
        }


        Self { size, 
            vector,
            min_temp, 
            max_temp, 
            temp,
        }

    }


    /// returns hotness based on max and min temp of fhr 
    pub fn hotness(&self, temp: ThermodynamicTemperature) -> f32 {

        let button_temp_degc = temp.get::<degree_celsius>();
        let min_temp_degc = self.min_temp.get::<degree_celsius>();
        let max_temp_degc = self.max_temp.get::<degree_celsius>();

        let hotness: f64 = 
            (button_temp_degc - min_temp_degc)/(max_temp_degc- min_temp_degc);

        return hotness as f32;
    }
    /// gets the size of the widget 
    pub fn size(&self) -> Vec2 {

        self.size.clone()
    }
}

impl Widget for SinglePipeColourBlueRedTempSensitive {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {

        let size = self.size();
        let (response, painter) = ui.allocate_painter(
            size, Sense::hover()
        );
        let pipe_hotness = 
            self.hotness(self.temp);

        let pipe_colour = hot_to_cold_colour_mark_1(
            pipe_hotness
        );
        // let colour = 
        let width = 20.0;

        let stroke = Stroke::new(width, pipe_colour);

        // get coordinates based on center
        let rect = response.rect;
        let pipe_ctr = rect.center();

        let delta_x = &self.vector.x;
        let delta_y = &self.vector.y;

        painter.line_segment(
            [pipe_ctr - vec2(0.50*delta_x, 0.50*delta_y), 
            pipe_ctr + vec2(0.50*delta_x, 0.50*delta_y)], 
            stroke
        );
        response
    }
}


pub struct SinglePipeColourBlackRedTempSensitive {
    size: Vec2,
    vector: Vec2,
    min_temp: ThermodynamicTemperature,
    max_temp: ThermodynamicTemperature,
    temp: ThermodynamicTemperature,
}

impl SinglePipeColourBlackRedTempSensitive {

    /// still need to correct for minimum size
    ///
    /// note that you will put in a vector from start to end 
    /// like how many pixels in x,y direction you want to go 
    /// then the pipe will autosize everything
    pub fn new(vector: Vec2,
        min_temp: ThermodynamicTemperature,
        max_temp: ThermodynamicTemperature,
        temp: ThermodynamicTemperature,) -> Self {

        let min_width = 20.0;

        // now the size here
        let mut size = vector;

        if size.x <= min_width {
            size.x = min_width
        };

        if size.y <= min_width {
            size.y = min_width
        }


        Self { size, 
            vector,
            min_temp, 
            max_temp, 
            temp,
        }

    }


    /// returns hotness based on max and min temp of fhr 
    pub fn hotness(&self, temp: ThermodynamicTemperature) -> f32 {

        let button_temp_degc = temp.get::<degree_celsius>();
        let min_temp_degc = self.min_temp.get::<degree_celsius>();
        let max_temp_degc = self.max_temp.get::<degree_celsius>();

        let hotness: f64 = 
            (button_temp_degc - min_temp_degc)/(max_temp_degc- min_temp_degc);

        return hotness as f32;
    }
    /// gets the size of the widget 
    pub fn size(&self) -> Vec2 {

        self.size.clone()
    }
}

impl Widget for SinglePipeColourBlackRedTempSensitive {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {

        let size = self.size();
        let (response, painter) = ui.allocate_painter(
            size, Sense::hover()
        );
        let pipe_hotness = 
            self.hotness(self.temp);

        let pipe_colour = hot_to_cold_colour_mark_2(
            pipe_hotness
        );
        // let colour = 
        let width = 20.0;

        let stroke = Stroke::new(width, pipe_colour);

        // get coordinates based on center
        let rect = response.rect;
        let pipe_ctr = rect.center();

        let delta_x = &self.vector.x;
        let delta_y = &self.vector.y;

        painter.line_segment(
            [pipe_ctr - vec2(0.50*delta_x, 0.50*delta_y), 
            pipe_ctr + vec2(0.50*delta_x, 0.50*delta_y)], 
            stroke
        );
        response
    }
}


pub struct SinglePipeColourBlueWhiteQualitySensitive {
    size: Vec2,
    vector: Vec2,
    quality: f64,
}

impl SinglePipeColourBlueWhiteQualitySensitive {

    /// still need to correct for minimum size
    ///
    /// note that you will put in a vector from start to end 
    /// like how many pixels in x,y direction you want to go 
    /// then the pipe will autosize everything
    pub fn new(vector: Vec2,
        quality: f64,) -> Self {

        let min_width = 20.0;

        // now the size here
        let mut size = vector;

        if size.x <= min_width {
            size.x = min_width
        };

        if size.y <= min_width {
            size.y = min_width
        }


        Self { size, 
            vector,
            quality,
        }

    }


    /// gets the size of the widget 
    pub fn size(&self) -> Vec2 {

        self.size.clone()
    }
}

impl Widget for SinglePipeColourBlueWhiteQualitySensitive {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {

        let size = self.size();
        let (response, painter) = ui.allocate_painter(
            size, Sense::hover()
        );
        let pipe_quality = self.quality;

        let pipe_colour = steam_quality_colour_mark_1(
            pipe_quality as f32
        );
        // let colour = 
        let width = 20.0;

        let stroke = Stroke::new(width, pipe_colour);

        // get coordinates based on center
        let rect = response.rect;
        let pipe_ctr = rect.center();

        let delta_x = &self.vector.x;
        let delta_y = &self.vector.y;

        painter.line_segment(
            [pipe_ctr - vec2(0.50*delta_x, 0.50*delta_y), 
            pipe_ctr + vec2(0.50*delta_x, 0.50*delta_y)], 
            stroke
        );
        response
    }
}


