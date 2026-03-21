use egui::{epaint::PathShape, vec2, Color32, Pos2, Sense, Stroke, Vec2, Widget};
use egui::epaint::CubicBezierShape;
use uom::si::ratio::ratio;
use uom::si::{f64::*, thermodynamic_temperature::degree_celsius};

use super::hot_to_cold_colour_mark_1;

#[derive(Debug,Clone, PartialEq)]
pub struct TurbineWidget {

    size: Vec2,

    theta: Angle,
}

// for the turbine widget, this is basically a collection of rectangles 
// with strokes in there. 
//
// These strokes will vary position with sine wave functions
impl TurbineWidget {

    pub fn new(size: Vec2, theta: Angle) -> Self {
        Self { size, theta  }
    }

    /// gets the size of the widget 
    pub fn size(&self) -> Vec2 {

        self.size.clone()
    }

    pub fn set_theta(&mut self, 
        omega: AngularVelocity,
        simulation_time: Time){
        self.theta = (omega * simulation_time).into();
    }

    pub fn get_theta(&self) -> Angle {
        self.theta
    }

}

impl Widget for TurbineWidget {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {
        let size = self.size();
        let (response, painter) = ui.allocate_painter(
            size, Sense::hover()
        );
        // this will be the main bulk of code

        let rect = response.rect;
        let c = rect.center();

        let rect_x = rect.width();
        let rect_y = rect.height();
        // we start with rectangles

        let turbine_blade_thickness = rect_x * 0.3;
        let turbine_radius_mid_blade = rect_y;


        let turbine_blade_colour = Color32::LIGHT_GRAY;


        let turbine_blade_stroke = Stroke::new(
            turbine_blade_thickness, 
            turbine_blade_colour
        );

        let turbine_center: Pos2 = c;

        painter.line_segment(
            [turbine_center - vec2(-0.20*turbine_blade_thickness, turbine_radius_mid_blade), 
            turbine_center + vec2(0.20*turbine_blade_thickness, turbine_radius_mid_blade)], 
            turbine_blade_stroke
        );

        // this is the spinning part
        let turbine_rotor_ratio = 0.20;
        let turbine_rotor_stroke = Stroke::new(
            turbine_blade_thickness * turbine_rotor_ratio, 
            Color32::DARK_GRAY
        );

        let rotor_mid_position_y: f64 = 
            rect_y as f64 * 
            self.get_theta().cos().get::<ratio>();

        let rotor_blade_center: Pos2 = 
            Pos2 { 
                x: turbine_center.x, 
                y: turbine_center.y + rotor_mid_position_y as f32            
            };




        // this makes a spinning blade
        //painter.line_segment(
        //    [turbine_center - vec2(-0.20*turbine_blade_thickness, stroke_mid_position_y as f32 * 1.1), 
        //    turbine_center + vec2(0.20*turbine_blade_thickness, stroke_mid_position_y as f32 * 0.9)], 
        //    turbine_rotor_stroke
        //);
        painter.line_segment(
            [rotor_blade_center - vec2(0.50*turbine_blade_thickness, 20.0), 
            rotor_blade_center + vec2(0.50*turbine_blade_thickness, 20.0)], 
            turbine_rotor_stroke
        );

        return response;
    }
}


