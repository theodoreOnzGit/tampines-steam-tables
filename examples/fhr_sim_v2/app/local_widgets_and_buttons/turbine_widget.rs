use egui::{epaint::PathShape, vec2, Color32, Pos2, Sense, Stroke, Vec2, Widget};
use egui::epaint::CubicBezierShape;
use uom::si::{f64::*, thermodynamic_temperature::degree_celsius};

use super::hot_to_cold_colour_mark_1;

#[derive(Debug,Clone, PartialEq)]
pub struct TurbineWidget {

    size: Vec2,
}

// for the turbine widget, this is basically a collection of rectangles 
// with strokes in there. 
//
// These strokes will vary position with sine wave functions
impl TurbineWidget {

    /// gets the size of the widget 
    pub fn size(&self) -> Vec2 {

        self.size.clone()
    }
}

impl Widget for TurbineWidget {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {
        let size = self.size();
        let (response, painter) = ui.allocate_painter(
            size, Sense::hover()
        );
        // this will be the main bulk of code




        return response;
    }
}
