use std::f64::consts::PI;

use egui::{epaint::PathShape, vec2, Color32, Pos2, Sense, Stroke, Vec2, Widget};
use egui::epaint::CubicBezierShape;
use uom::ConstZero;
use uom::si::angle::radian;
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

        let turbine_blade_axial_thickness = 0.1 * rect_x;
        let turbine_max_radius = 0.5*rect_y;


        let turbine_blade_colour = Color32::GRAY;


        let turbine_blade_stroke = Stroke::new(
            turbine_blade_axial_thickness, 
            turbine_blade_colour
        );

        let turbine_center: Pos2 = c;


        // this is the spinning part
        //
        // the rotor ratio is the relative ratio of the blade to the 
        // axial width of the turbine
        let turbine_rotor_ratio = 0.20;
        let turbine_rotor_stroke = Stroke::new(
            turbine_blade_axial_thickness * turbine_rotor_ratio, 
            Color32::BLACK
        );
        let turbine_rotor_stroke_white = Stroke::new(
            turbine_blade_axial_thickness * turbine_rotor_ratio, 
            Color32::WHITE
        );

        let turbine_num_blades = 20;
        let turbine_num_axial_sides = 5;

        let turbine_min_radius_pixels: f32 = 
            turbine_max_radius / turbine_num_axial_sides as f32;

        
        let paint_turbine_blade_set = |set_number: isize|{

            let turbine_radius: f32 = 
                turbine_min_radius_pixels 
                * (set_number - 1) as f32;


            let offset_factor_to_the_right = set_number as f32;

            // first the stroke for the turbine
            painter.line_segment(
                [turbine_center - vec2(
                    -turbine_blade_axial_thickness * offset_factor_to_the_right, 
                    turbine_radius
                ), 
                turbine_center + vec2(
                    turbine_blade_axial_thickness * offset_factor_to_the_right, 
                    turbine_radius
                )], 
                turbine_blade_stroke
            );
            for i in 0..turbine_num_blades {

                let theta_plus_phase_shift = 
                    self.get_theta()
                    + 
                    Angle::new::<radian>(i as f64 * (2.0 * PI)/(turbine_num_blades as f64))
                    ;

                let rotor_mid_position_y: f64 = 
                    0.9 * turbine_radius as f64 * 
                    theta_plus_phase_shift.cos().get::<ratio>();

                let rotor_blade_center: Pos2 = 
                    Pos2 { 
                        x: turbine_center.x + offset_factor_to_the_right * turbine_blade_axial_thickness, 
                        y: turbine_center.y + rotor_mid_position_y as f32            
                    };




                // this makes a spinning blade
                //painter.line_segment(
                //    [turbine_center - vec2(-0.20*turbine_blade_thickness, stroke_mid_position_y as f32 * 1.1), 
                //    turbine_center + vec2(0.20*turbine_blade_thickness, stroke_mid_position_y as f32 * 0.9)], 
                //    turbine_rotor_stroke

                //);
                // only paint moving blades if turbine 
                if theta_plus_phase_shift.sin() > Ratio::ZERO {



                    painter.line_segment(
                        [rotor_blade_center - vec2(
                            0.5*turbine_blade_axial_thickness, 
                            0.1 * turbine_radius
                        ), 
                        rotor_blade_center + vec2(
                            0.5*turbine_blade_axial_thickness, 
                            0.1 * turbine_radius
                        )], 
                        turbine_rotor_stroke
                    );

                    if i == 10 {
                        painter.line_segment(
                            [rotor_blade_center - vec2(
                                0.5*turbine_blade_axial_thickness, 
                                0.1 * turbine_radius
                            ), 
                            rotor_blade_center + vec2(
                                0.5*turbine_blade_axial_thickness, 
                                0.1 * turbine_radius
                            )], 
                            turbine_rotor_stroke_white
                        );

                    }
                }


                // end of blade painting
            }
        };

        //paint_turbine_blade_set(-1);
        paint_turbine_blade_set(-2);
        //for i in -5..=5 {
        //    paint_turbine_blade_set(i);
        //}


        //for i in 0..turbine_num_blades {

        //    let theta_plus_phase_shift = 
        //        self.get_theta()
        //        + 
        //        Angle::new::<radian>(i as f64 * (2.0 * PI)/(turbine_num_blades as f64))
        //        ;

        //    let rotor_mid_position_y: f64 = 
        //        0.4 * turbine_radius_mid_blade as f64 * 
        //        theta_plus_phase_shift.cos().get::<ratio>();

        //    let rotor_blade_center: Pos2 = 
        //        Pos2 { 
        //            x: turbine_center.x, 
        //            y: turbine_center.y + rotor_mid_position_y as f32            
        //        };




        //    // this makes a spinning blade
        //    //painter.line_segment(
        //    //    [turbine_center - vec2(-0.20*turbine_blade_thickness, stroke_mid_position_y as f32 * 1.1), 
        //    //    turbine_center + vec2(0.20*turbine_blade_thickness, stroke_mid_position_y as f32 * 0.9)], 
        //    //    turbine_rotor_stroke

        //    //);
        //    // only paint moving blades if turbine 
        //    if theta_plus_phase_shift.sin() > Ratio::ZERO {
        //        painter.line_segment(
        //            [rotor_blade_center - vec2(0.50*turbine_blade_thickness, 10.0), 
        //            rotor_blade_center + vec2(0.50*turbine_blade_thickness, 10.0)], 
        //            turbine_rotor_stroke
        //        );

        //        if i == 10 {
        //        painter.line_segment(
        //            [rotor_blade_center - vec2(0.50*turbine_blade_thickness, 10.0), 
        //            rotor_blade_center + vec2(0.50*turbine_blade_thickness, 10.0)], 
        //            turbine_rotor_stroke_red
        //        );

        //        }
        //    }
        //}

        return response;
    }
}


