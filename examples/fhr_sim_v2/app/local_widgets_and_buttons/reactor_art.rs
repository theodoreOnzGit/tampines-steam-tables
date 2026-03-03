
use egui::{epaint::{CubicBezierShape, PathShape}, vec2, Color32, Pos2, Rect, Sense, Stroke, Ui, Vec2};

use super::hot_to_cold_colour_mark_1;

// first thing, reactor
//
// it needs a top left coordinate first
// these will make the rectangle
//
// i think for easy gui connections, I can start making classes 
// and stuff for easy connection
pub fn fhr_reactor_vessel_prototype(ui: &mut Ui,
    user_rect: egui::Rect,
    mut left_control_rod_insertion_frac: f32,
    mut right_control_rod_insertion_frac: f32,
    ){

    // make a new painter first 
    //

    if left_control_rod_insertion_frac > 1.0 {
        left_control_rod_insertion_frac = 1.0;
    } else if left_control_rod_insertion_frac < 0.0 {
        left_control_rod_insertion_frac = 0.0;
    };
    if right_control_rod_insertion_frac > 1.0 {
        right_control_rod_insertion_frac = 1.0;
    } else if right_control_rod_insertion_frac < 0.0 {
        right_control_rod_insertion_frac = 0.0;
    };



    // top_left
    let left_most_side = user_rect.left();
    let top_most_side = user_rect.top();
    let right_most_side = user_rect.right();
    let bottom_most_side = user_rect.bottom();

    let max_height_y = (top_most_side - bottom_most_side).abs();
    let max_width_x = (left_most_side - right_most_side).abs();

    let ui_rectangle: Rect = ui.min_rect();
    let breadth = ui_rectangle.right();
    let height = ui_rectangle.bottom();

    // the size here is the size of the painter
    let size = Vec2::new(breadth, height);

    let (response, painter) = ui.allocate_painter(
        size, Sense::hover()
    );
    let _response_rect = response.rect;
    
    // what I want to do now is to shift the 
    // response rectangle

    let _shift: Pos2 = user_rect.min;

    // the rect here is 
    let rect: egui::Rect = 
        egui:: Rect {
            min: user_rect.min,
            max: user_rect.max,
        };
    
    let c = rect.center();
    // circle radius is r
    let _r = rect.width() / 2.0 - 1.0;
    let color = Color32::from_gray(128);
    let stroke = Stroke::new(1.0, color);

    let reactor_half_length_y = max_height_y * 0.5 * 0.8;
    let reactor_half_width_x = max_width_x * 0.5 * 0.8;
    //painter.circle_stroke(c, r, stroke);
    
    // now let's paint the reactor first
    //painter.line_segment(
    //    [c + vec2(reactor_half_width_x, reactor_half_length_y), 
    //    c + vec2(reactor_half_width_x, -reactor_half_length_y)], 
    //    stroke
    //);
    //painter.line_segment(
    //    [c + vec2(-reactor_half_width_x, reactor_half_length_y), 
    //    c + vec2(-reactor_half_width_x, -reactor_half_length_y)], 
    //    stroke
    //);
    //painter.line_segment([c, c + r * Vec2::angled(TAU * 1.0 / 8.0)], stroke);
    //painter.line_segment([c, c + r * Vec2::angled(TAU * 3.0 / 8.0)], stroke);
    //
    // if wanting filled shapes:
    // https://docs.rs/epaint/0.31.1/epaint/struct.PathShape.html
    //
    // also Beizier curve
    // https://docs.rs/epaint/0.31.1/epaint/struct.QuadraticBezierShape.html
    // https://docs.rs/epaint/0.31.1/epaint/struct.CubicBezierShape.html
    //
    // with these filled curves, I can do a shape which 
    // for which I can then use the add method to add shapes 
    // to the painter:
    // 
    // https://docs.rs/egui/latest/egui/struct.Painter.html
    //
    // so to get the shape, i have a complex polygon
    // The points are

    

    // bottom inlet
    let fhr_coolant_inlet_bottom_left = c + vec2(-0.10*reactor_half_width_x, reactor_half_length_y* 0.9);
    let fhr_coolant_inlet_bottom_right = c + vec2(0.10*reactor_half_width_x, reactor_half_length_y* 0.9);

    // core part
    let fhr_core_inlet_bottom_left = c + vec2(-0.10*reactor_half_width_x, reactor_half_length_y* 0.65);
    let fhr_core_inlet_bottom_right = c + vec2(0.10*reactor_half_width_x, reactor_half_length_y* 0.65);
    let fhr_core_fat_bottom_left = c + vec2(-0.50*reactor_half_width_x, reactor_half_length_y* 0.45);
    let fhr_core_fat_bottom_right = c + vec2(0.50*reactor_half_width_x, reactor_half_length_y* 0.45);

    let fhr_core_fat_top_left = c + vec2(-0.50*reactor_half_width_x, -reactor_half_length_y* 0.45);
    let fhr_core_fat_top_right = c + vec2(0.50*reactor_half_width_x, -reactor_half_length_y* 0.45);
    let fhr_core_outlet_top_left = c + vec2(-0.10*reactor_half_width_x, -reactor_half_length_y* 0.65);
    let fhr_core_outlet_top_right = c + vec2(0.10*reactor_half_width_x, -reactor_half_length_y* 0.65);

    // top outlet
    let fhr_coolant_outlet_top_left = c + vec2(-0.10*reactor_half_width_x, -reactor_half_length_y* 0.9);
    let fhr_coolant_outlet_top_right = c + vec2(0.10*reactor_half_width_x, -reactor_half_length_y* 0.9);

    // colour fill
    let hotness: f32 = 0.1;
    let coolant_fill = hot_to_cold_colour_mark_1(hotness);

    // draw clockwise
    let core_bottom_points = vec![
        //fhr_coolant_inlet_bottom_left,
        fhr_core_inlet_bottom_left,
        fhr_core_fat_bottom_left,
        //fhr_core_fat_top_left,
        //fhr_core_outlet_top_left,
        //fhr_coolant_outlet_top_left,
        //fhr_coolant_outlet_top_right,
        //fhr_core_outlet_top_right,
        //fhr_core_fat_top_right,
        fhr_core_fat_bottom_right,
        fhr_core_inlet_bottom_right,
        //fhr_coolant_inlet_bottom_right,
    ];
    let core_bottom_inlet_points = vec![
        fhr_coolant_inlet_bottom_left,
        fhr_core_inlet_bottom_left,
        //fhr_core_fat_bottom_left,
        //fhr_core_fat_top_left,
        //fhr_core_outlet_top_left,
        //fhr_coolant_outlet_top_left,
        //fhr_coolant_outlet_top_right,
        //fhr_core_outlet_top_right,
        //fhr_core_fat_top_right,
        //fhr_core_fat_bottom_right,
        fhr_core_inlet_bottom_right,
        fhr_coolant_inlet_bottom_right,
    ];
    let core_mid_points = vec![
        //fhr_coolant_inlet_bottom_left,
        //fhr_core_inlet_bottom_left,
        fhr_core_fat_bottom_left,
        fhr_core_fat_top_left,
        //fhr_core_outlet_top_left,
        //fhr_coolant_outlet_top_left,
        //fhr_coolant_outlet_top_right,
        //fhr_core_outlet_top_right,
        fhr_core_fat_top_right,
        fhr_core_fat_bottom_right,
        //fhr_core_inlet_bottom_right,
        //fhr_coolant_inlet_bottom_right,
    ];
    let core_top_points = vec![
        //fhr_coolant_inlet_bottom_left,
        //fhr_core_inlet_bottom_left,
        //fhr_core_fat_bottom_left,
        fhr_core_fat_top_left,
        fhr_core_outlet_top_left,
        //fhr_coolant_outlet_top_left,
        //fhr_coolant_outlet_top_right,
        fhr_core_outlet_top_right,
        fhr_core_fat_top_right,
        //fhr_core_fat_bottom_right,
        //fhr_core_inlet_bottom_right,
        //fhr_coolant_inlet_bottom_right,
    ];
    let core_outlet_points = vec![
        //fhr_coolant_inlet_bottom_left,
        //fhr_core_inlet_bottom_left,
        //fhr_core_fat_bottom_left,
        //fhr_core_fat_top_left,
        fhr_core_outlet_top_left,
        fhr_coolant_outlet_top_left,
        fhr_coolant_outlet_top_right,
        fhr_core_outlet_top_right,
        //fhr_core_fat_top_right,
        //fhr_core_fat_bottom_right,
        //fhr_core_inlet_bottom_right,
        //fhr_coolant_inlet_bottom_right,
    ];

    // fhr metal container (grey colour)
    //
    // we will use a cubic Beizier curve


    let reactor_box_top_left = 
        c + vec2(-reactor_half_width_x, -reactor_half_length_y);
    let reactor_box_bottom_left = 
        c + vec2(-reactor_half_width_x, reactor_half_length_y);
    let reactor_box_top_right = 
        c + vec2(reactor_half_width_x, -reactor_half_length_y);
    let reactor_box_bottom_right = 
        c + vec2(reactor_half_width_x, reactor_half_length_y);

    
    let reactor_curved_edge_fraction = 0.55;

    let reactor_curved_edge_top_left = 
        c + vec2(-reactor_half_width_x, -reactor_curved_edge_fraction * reactor_half_length_y);
    let reactor_curved_edge_bottom_left = 
        c + vec2(-reactor_half_width_x, reactor_curved_edge_fraction * reactor_half_length_y);
    let reactor_curved_edge_top_right = 
        c + vec2(reactor_half_width_x, -reactor_curved_edge_fraction * reactor_half_length_y);
    let reactor_curved_edge_bottom_right = 
        c + vec2(reactor_half_width_x, reactor_curved_edge_fraction * reactor_half_length_y);

    let metal_fill = Color32::GRAY;

    let fhr_bottom_metal_pts = 
        [
        reactor_curved_edge_bottom_left,
        reactor_box_bottom_left,
        reactor_box_bottom_right,
        reactor_curved_edge_bottom_right
        ];
    let fhr_top_metal_pts = 
        [
        reactor_curved_edge_top_left,
        reactor_box_top_left,
        reactor_box_top_right,
        reactor_curved_edge_top_right
        ];

    let fhr_mid_metal_pts = 
        [
        reactor_curved_edge_bottom_left,
        reactor_curved_edge_top_left,
        reactor_curved_edge_top_right,
        reactor_curved_edge_bottom_right
        ];

    let fhr_bottom_metal_semicircle = 
        CubicBezierShape::from_points_stroke(fhr_bottom_metal_pts, 
            true, 
            metal_fill, 
            stroke);

    let fhr_top_metal_semicircle = 
        CubicBezierShape::from_points_stroke(fhr_top_metal_pts, 
            true, 
            metal_fill, 
            stroke);
    let fhr_mid_metal_rect = 
        PathShape::convex_polygon(fhr_mid_metal_pts.into(), 
            metal_fill, 
            stroke);

    // fhr metal vessel
    painter.add(fhr_bottom_metal_semicircle);
    painter.add(fhr_top_metal_semicircle);
    painter.add(fhr_mid_metal_rect);


    // inner graphite reflector 

    let graphite_width_fraction = 0.8;

    let reflector_box_top_left = 
        c + vec2(-reactor_half_width_x * graphite_width_fraction, -reactor_half_length_y * graphite_width_fraction);
    let reflector_box_bottom_left = 
        c + vec2(-reactor_half_width_x * graphite_width_fraction, reactor_half_length_y * graphite_width_fraction);
    let reflector_box_top_right = 
        c + vec2(reactor_half_width_x * graphite_width_fraction, -reactor_half_length_y * graphite_width_fraction);
    let reflector_box_bottom_right = 
        c + vec2(reactor_half_width_x * graphite_width_fraction, reactor_half_length_y * graphite_width_fraction);

    
    let reflector_curved_edge_fraction = 0.55;

    let reflector_curved_edge_top_left = 
        c + vec2(-reactor_half_width_x * graphite_width_fraction, -reflector_curved_edge_fraction * reactor_half_length_y);
    let reflector_curved_edge_bottom_left = 
        c + vec2(-reactor_half_width_x * graphite_width_fraction, reflector_curved_edge_fraction * reactor_half_length_y);
    let reflector_curved_edge_top_right = 
        c + vec2(reactor_half_width_x * graphite_width_fraction, -reflector_curved_edge_fraction * reactor_half_length_y);
    let reflector_curved_edge_bottom_right = 
        c + vec2(reactor_half_width_x * graphite_width_fraction, reflector_curved_edge_fraction * reactor_half_length_y);


    let reflector_bottom_graphite_pts = 
        [
        reflector_curved_edge_bottom_left,
        reflector_box_bottom_left,
        reflector_box_bottom_right,
        reflector_curved_edge_bottom_right
        ];
    let reflector_top_graphite_pts = 
        [
        reflector_curved_edge_top_left,
        reflector_box_top_left,
        reflector_box_top_right,
        reflector_curved_edge_top_right
        ];

    let reflector_mid_graphite_pts = 
        [
        reflector_curved_edge_bottom_left,
        reflector_curved_edge_top_left,
        reflector_curved_edge_top_right,
        reflector_curved_edge_bottom_right
        ];
    let graphite_fill = Color32::BLACK;

    let graphite_stroke = Stroke::new(1.0, graphite_fill);
    let reflector_bottom_graphite_semicircle = 
        CubicBezierShape::from_points_stroke(reflector_bottom_graphite_pts, 
            true, 
            graphite_fill, 
            graphite_stroke);

    let reflector_top_graphite_semicircle = 
        CubicBezierShape::from_points_stroke(reflector_top_graphite_pts, 
            true, 
            graphite_fill, 
            graphite_stroke);
    let reflector_mid_graphite_rect = 
        PathShape::convex_polygon(reflector_mid_graphite_pts.into(), 
            graphite_fill, 
            graphite_stroke);

    // fhr reflector graphite
    painter.add(reflector_bottom_graphite_semicircle);
    painter.add(reflector_top_graphite_semicircle);
    painter.add(reflector_mid_graphite_rect);
    
    let coolant_stroke = Stroke::new(1.0, coolant_fill);
    // fhr coolant 
    let fhr_core_bottom_coolant_shape = 
        PathShape::convex_polygon(
            core_bottom_points, 
            coolant_fill, 
            coolant_stroke);
    let fhr_core_inlet_coolant_shape = 
        PathShape::convex_polygon(
            core_bottom_inlet_points, 
            coolant_fill, 
            coolant_stroke);
    let fhr_core_mid_coolant_shape = 
        PathShape::convex_polygon(
            core_mid_points, 
            coolant_fill, 
            stroke);
    let fhr_core_top_coolant_shape = 
        PathShape::convex_polygon(
            core_top_points, 
            coolant_fill, 
            coolant_stroke);
    let fhr_core_outlet_coolant_shape = 
        PathShape::convex_polygon(
            core_outlet_points, 
            coolant_fill, 
            coolant_stroke);



    painter.add(fhr_core_bottom_coolant_shape);
    painter.add(fhr_core_inlet_coolant_shape);
    painter.add(fhr_core_top_coolant_shape);
    painter.add(fhr_core_outlet_coolant_shape);
    painter.add(fhr_core_mid_coolant_shape);

    // now for pebble bed 
    //
    let fhr_width = max_width_x * 0.4;
    let fhr_height = max_height_y * 0.8;
    let pebble_radius = fhr_width * 0.08;
    let core_radius = pebble_radius * 0.8;
    let pebble_ctr = c;

    painter.circle_filled(pebble_ctr, pebble_radius, Color32::BLACK);
    painter.circle_filled(pebble_ctr, core_radius, Color32::DARK_RED);

    let pebble_centers = vec![
        c + vec2(2.0*pebble_radius,0.1*pebble_radius),
        c + vec2(1.0*pebble_radius,-0.5*pebble_radius),
        c + vec2(4.0*pebble_radius,-0.3*pebble_radius),
        c + vec2(-2.0*pebble_radius,0.1*pebble_radius),
        c + vec2(-1.0*pebble_radius,-0.5*pebble_radius),
        c + vec2(-4.0*pebble_radius,-0.3*pebble_radius),
        c + vec2(2.0*pebble_radius,1.1*pebble_radius),
        c + vec2(1.0*pebble_radius,-1.5*pebble_radius),
        c + vec2(4.0*pebble_radius,-2.3*pebble_radius),
        c + vec2(-2.0*pebble_radius,1.1*pebble_radius),
        c + vec2(-1.0*pebble_radius,-2.5*pebble_radius),
        c + vec2(-4.0*pebble_radius,-2.3*pebble_radius),
        c + vec2(2.0*pebble_radius,1.1*pebble_radius),
        c + vec2(3.0*pebble_radius,-1.5*pebble_radius),
        c + vec2(5.0*pebble_radius,-2.3*pebble_radius),
        c + vec2(-2.0*pebble_radius,1.1*pebble_radius),
        c + vec2(-3.0*pebble_radius,-2.5*pebble_radius),
        c + vec2(-5.0*pebble_radius,-2.3*pebble_radius),
        c + vec2(-5.0*pebble_radius,2.3*pebble_radius),
        c + vec2(5.0*pebble_radius,2.3*pebble_radius),
        c + vec2(-4.2*pebble_radius,1.3*pebble_radius),
        c + vec2(4.0*pebble_radius,1.4*pebble_radius),
        c + vec2(-0.2*pebble_radius,1.3*pebble_radius),
        c + vec2(0.0*pebble_radius,1.8*pebble_radius),
    ];

    // add another list but transpose upwards 
    let mut pebble_centres_bottom: Vec<Pos2> = pebble_centers.clone();
    let mut pebble_centres_top: Vec<Pos2> = pebble_centers.clone();

    for (i,pebble_center) in pebble_centers.iter().enumerate() {
        pebble_centres_bottom[i] = 
            *pebble_center + vec2(0.0, fhr_height * 0.1);
    }
    for (i,pebble_center) in pebble_centers.iter().enumerate() {
        pebble_centres_top[i] = 
            *pebble_center + vec2(0.0, -fhr_height * 0.1);
    }


    for pebble_center in pebble_centers.iter(){
        painter.circle_filled(*pebble_center, pebble_radius, Color32::BLACK);
        painter.circle_filled(*pebble_center, core_radius, Color32::DARK_RED);
    }
    
    for pebble_center in pebble_centres_bottom.iter(){
        painter.circle_filled(*pebble_center, pebble_radius, Color32::BLACK);
        painter.circle_filled(*pebble_center, core_radius, Color32::DARK_RED);
    }
    for pebble_center in pebble_centres_top.iter(){
        painter.circle_filled(*pebble_center, pebble_radius, Color32::BLACK);
        painter.circle_filled(*pebble_center, core_radius, Color32::DARK_RED);
    }

    // next, downcomers

    // left downcomer inlet
    let left_downcomer_inlet_bottom_pt = 
        fhr_coolant_inlet_bottom_left 
        + vec2(0.0, -reactor_half_length_y * 0.04);

    let left_downcomer_inlet_top_pt = 
        fhr_coolant_inlet_bottom_left 
        + vec2(0.0, -reactor_half_length_y * 0.12);

    let left_downcomer_inlet_mid_bottom_pt = 
        fhr_coolant_inlet_bottom_left 
        + vec2(-reactor_half_width_x * 0.65, -reactor_half_length_y * 0.16);

    let left_downcomer_inlet_mid_top_pt = 
        fhr_coolant_inlet_bottom_left 
        + vec2(-reactor_half_width_x * 0.6, -reactor_half_length_y * 0.22);

    // left downcomer mid rectangle
    //
    let left_downcomer_mid_rect_bottom_left =
        reactor_curved_edge_bottom_left 
        + vec2(reactor_half_width_x *0.06, 0.0);

    let left_downcomer_mid_rect_bottom_right =
        reactor_curved_edge_bottom_left 
        + vec2(reactor_half_width_x *0.16, 0.0);


    let left_downcomer_mid_rect_top_left =
        reactor_curved_edge_top_left 
        + vec2(reactor_half_width_x *0.06, 0.0);

    let left_downcomer_mid_rect_top_right =
        reactor_curved_edge_top_left 
        + vec2(reactor_half_width_x *0.16, 0.0);

    // left downcomer outlet

    let left_downcomer_outlet_top_pt = 
        fhr_coolant_outlet_top_left 
        + vec2(0.0, reactor_half_length_y * 0.04);

    let left_downcomer_outlet_bottom_pt = 
        fhr_coolant_outlet_top_left 
        + vec2(0.0, reactor_half_length_y * 0.12);

    let left_downcomer_outlet_mid_bottom_pt = 
        fhr_coolant_outlet_top_left 
        + vec2(-reactor_half_width_x * 0.6, reactor_half_length_y * 0.22);

    let left_downcomer_outlet_mid_top_pt = 
        fhr_coolant_outlet_top_left 
        + vec2(-reactor_half_width_x * 0.65, reactor_half_length_y * 0.16);



    let downcomer_inlet_left_1_pts = 
        vec![
        left_downcomer_inlet_bottom_pt,
        left_downcomer_inlet_mid_bottom_pt,
        left_downcomer_inlet_mid_top_pt,
        left_downcomer_inlet_top_pt
        ];

    let downcomer_inlet_left_2_pts = 
        vec![
        left_downcomer_mid_rect_bottom_left,
        left_downcomer_mid_rect_bottom_right,
        left_downcomer_inlet_mid_top_pt,
        left_downcomer_inlet_mid_bottom_pt,
        ];


    let downcomer_left_mid_pts = 
        vec![
        left_downcomer_mid_rect_bottom_left,
        left_downcomer_mid_rect_top_left,
        left_downcomer_mid_rect_top_right,
        left_downcomer_mid_rect_bottom_right,
        ];

    let downcomer_outlet_left_1_pts = 
        vec![
        left_downcomer_outlet_bottom_pt,
        left_downcomer_outlet_mid_bottom_pt,
        left_downcomer_outlet_mid_top_pt,
        left_downcomer_outlet_top_pt
        ];

    let downcomer_outlet_left_2_pts = 
        vec![
        left_downcomer_mid_rect_top_left,
        left_downcomer_outlet_mid_top_pt,
        left_downcomer_outlet_mid_bottom_pt,
        left_downcomer_mid_rect_top_right,
        ];

    let left_downcomer_inlet_1_shape = 
        PathShape::convex_polygon(
            downcomer_inlet_left_1_pts, 
            coolant_fill, 
            coolant_stroke);
    let left_downcomer_inlet_2_shape = 
        PathShape::convex_polygon(
            downcomer_inlet_left_2_pts, 
            coolant_fill, 
            coolant_stroke);
    let left_downcomer_mid_shape = 
        PathShape::convex_polygon(
            downcomer_left_mid_pts, 
            coolant_fill, 
            coolant_stroke);
    let left_downcomer_outlet_1_shape = 
        PathShape::convex_polygon(
            downcomer_outlet_left_1_pts, 
            coolant_fill, 
            coolant_stroke);
    let left_downcomer_outlet_2_shape = 
        PathShape::convex_polygon(
            downcomer_outlet_left_2_pts, 
            coolant_fill, 
            coolant_stroke);

    painter.add(left_downcomer_inlet_1_shape);
    painter.add(left_downcomer_inlet_2_shape);
    painter.add(left_downcomer_mid_shape);
    painter.add(left_downcomer_outlet_1_shape);
    painter.add(left_downcomer_outlet_2_shape);

    // right downcomer 

    // right downcomer inlet
    let right_downcomer_inlet_bottom_pt = 
        fhr_coolant_inlet_bottom_right 
        + vec2(0.0, -reactor_half_length_y * 0.04);

    let right_downcomer_inlet_top_pt = 
        fhr_coolant_inlet_bottom_right 
        + vec2(0.0, -reactor_half_length_y * 0.12);

    let right_downcomer_inlet_mid_bottom_pt = 
        fhr_coolant_inlet_bottom_right 
        + vec2(reactor_half_width_x * 0.65, -reactor_half_length_y * 0.16);

    let right_downcomer_inlet_mid_top_pt = 
        fhr_coolant_inlet_bottom_right 
        + vec2(reactor_half_width_x * 0.6, -reactor_half_length_y * 0.22);

    // right downcomer mid rectangle
    //
    let right_downcomer_mid_rect_bottom_left =
        reactor_curved_edge_bottom_right 
        + vec2(-reactor_half_width_x *0.16, 0.0);

    let right_downcomer_mid_rect_bottom_right =
        reactor_curved_edge_bottom_right 
        + vec2(-reactor_half_width_x *0.06, 0.0);


    let right_downcomer_mid_rect_top_left =
        reactor_curved_edge_top_right 
        + vec2(-reactor_half_width_x *0.16, 0.0);

    let right_downcomer_mid_rect_top_right =
        reactor_curved_edge_top_right 
        + vec2(-reactor_half_width_x *0.06, 0.0);

    // right downcomer outlet

    let right_downcomer_outlet_top_pt = 
        fhr_coolant_outlet_top_right 
        + vec2(0.0, reactor_half_length_y * 0.04);

    let right_downcomer_outlet_bottom_pt = 
        fhr_coolant_outlet_top_right 
        + vec2(0.0, reactor_half_length_y * 0.12);

    let right_downcomer_outlet_mid_bottom_pt = 
        fhr_coolant_outlet_top_right 
        + vec2(reactor_half_width_x * 0.6, reactor_half_length_y * 0.22);

    let right_downcomer_outlet_mid_top_pt = 
        fhr_coolant_outlet_top_right 
        + vec2(reactor_half_width_x * 0.65, reactor_half_length_y * 0.16);

    let downcomer_inlet_right_1_pts = 
        vec![
        right_downcomer_inlet_bottom_pt,
        right_downcomer_inlet_mid_bottom_pt,
        right_downcomer_inlet_mid_top_pt,
        right_downcomer_inlet_top_pt
        ];

    let downcomer_inlet_right_2_pts = 
        vec![
        right_downcomer_inlet_mid_top_pt,
        right_downcomer_mid_rect_bottom_left,
        right_downcomer_mid_rect_bottom_right,
        right_downcomer_inlet_mid_bottom_pt,
        ];


    let downcomer_right_mid_pts = 
        vec![
        right_downcomer_mid_rect_bottom_left,
        right_downcomer_mid_rect_top_left,
        right_downcomer_mid_rect_top_right,
        right_downcomer_mid_rect_bottom_right,
        ];

    let downcomer_outlet_right_1_pts = 
        vec![
        right_downcomer_outlet_bottom_pt,
        right_downcomer_outlet_mid_bottom_pt,
        right_downcomer_outlet_mid_top_pt,
        right_downcomer_outlet_top_pt
        ];

    let downcomer_outlet_right_2_pts = 
        vec![
        right_downcomer_outlet_mid_top_pt,
        right_downcomer_mid_rect_top_right,
        right_downcomer_mid_rect_top_left,
        right_downcomer_outlet_mid_bottom_pt,
        ];

    let right_downcomer_inlet_1_shape = 
        PathShape::convex_polygon(
            downcomer_inlet_right_1_pts, 
            coolant_fill, 
            coolant_stroke);
    let right_downcomer_inlet_2_shape = 
        PathShape::convex_polygon(
            downcomer_inlet_right_2_pts, 
            coolant_fill, 
            coolant_stroke);
    let right_downcomer_mid_shape = 
        PathShape::convex_polygon(
            downcomer_right_mid_pts, 
            coolant_fill, 
            coolant_stroke);
    let right_downcomer_outlet_1_shape = 
        PathShape::convex_polygon(
            downcomer_outlet_right_1_pts, 
            coolant_fill, 
            coolant_stroke);
    let right_downcomer_outlet_2_shape = 
        PathShape::convex_polygon(
            downcomer_outlet_right_2_pts, 
            coolant_fill, 
            coolant_stroke);

    painter.add(right_downcomer_inlet_1_shape);
    painter.add(right_downcomer_inlet_2_shape);
    painter.add(right_downcomer_mid_shape);
    painter.add(right_downcomer_outlet_1_shape);
    painter.add(right_downcomer_outlet_2_shape);

    // now control rods

    let cr_channel_length_ratio = 1.0;
    let cr_channel_width_ratio = 0.08;

    let cr_left_ref_pt = 
        reflector_curved_edge_top_left 
        + vec2(reactor_half_width_x * 0.15, 0.0);

    let left_cr_channel_top_left = 
        cr_left_ref_pt 
        + vec2(-reactor_half_width_x * cr_channel_width_ratio, 0.0);

    let left_cr_channel_top_right = 
        cr_left_ref_pt 
        + vec2(reactor_half_width_x * cr_channel_width_ratio, 0.0);

    let left_cr_channel_bottom_left = 
        cr_left_ref_pt 
        + vec2(-reactor_half_width_x * 0.08, 
            reactor_half_length_y * cr_channel_length_ratio);

    let left_cr_channel_bottom_right = 
        cr_left_ref_pt 
        + vec2(reactor_half_width_x * 0.08, 
            reactor_half_length_y * cr_channel_length_ratio);

    let cr_channel_fill = Color32::LIGHT_BLUE;




    let left_cr_channel_pts = 
        vec![
        left_cr_channel_top_left,
        left_cr_channel_top_right,
        left_cr_channel_bottom_right,
        left_cr_channel_bottom_left,
        ];

    let cr_left_channel_shape = 
        PathShape::convex_polygon(
            left_cr_channel_pts, 
            cr_channel_fill, 
            coolant_stroke);


    let cr_right_ref_pt = 
        reflector_curved_edge_top_right 
        + vec2(-reactor_half_width_x * 0.15, 0.0);

    let right_cr_channel_top_left = 
        cr_right_ref_pt 
        + vec2(-reactor_half_width_x * cr_channel_width_ratio, 0.0);

    let right_cr_channel_top_right = 
        cr_right_ref_pt 
        + vec2(reactor_half_width_x * cr_channel_width_ratio, 0.0);

    let right_cr_channel_bottom_left = 
        cr_right_ref_pt 
        + vec2(-reactor_half_width_x * 0.08, 
            reactor_half_length_y * cr_channel_length_ratio);

    let right_cr_channel_bottom_right = 
        cr_right_ref_pt 
        + vec2(reactor_half_width_x * 0.08, 
            reactor_half_length_y * cr_channel_length_ratio);
    let right_cr_channel_pts = 
        vec![
        right_cr_channel_top_left,
        right_cr_channel_top_right,
        right_cr_channel_bottom_right,
        right_cr_channel_bottom_left,
        ];

    let cr_right_channel_shape = 
        PathShape::convex_polygon(
            right_cr_channel_pts, 
            cr_channel_fill, 
            coolant_stroke);

    painter.add(cr_left_channel_shape);
    painter.add(cr_right_channel_shape);

    let cr_width_ratio = 0.08;
    let cr_colour = Color32::DARK_GRAY;


    let cr_rod_stroke = Stroke::new(
        cr_width_ratio * reactor_half_width_x, 
        cr_colour
    );

    let cr_length = reactor_half_length_y * 0.88;

    let cr_left_centre = cr_left_ref_pt 
        + vec2(0.0, cr_length*left_control_rod_insertion_frac - cr_length*0.9);

    let cr_right_centre = cr_right_ref_pt 
        + vec2(0.0, cr_length*right_control_rod_insertion_frac - cr_length*0.9);


    painter.line_segment(
        [cr_left_centre - vec2(0.0, cr_length), 
        cr_left_centre + vec2(0.0, cr_length)], 
        cr_rod_stroke
    );
    painter.line_segment(
        [cr_right_centre - vec2(0.0, cr_length), 
        cr_right_centre + vec2(0.0, cr_length)], 
        cr_rod_stroke
    );

    //let label = egui::Label::new("dummy label");

    //let label_top_left = reactor_box_bottom_left;
    //let label_top_right = reactor_box_top_right;
    //let label_bottom_right = label_top_right 
    //    + vec2(0.0, -reactor_half_length_y*0.08);
    //let label_rect = egui::Rect {
    //    min: label_top_left,
    //    max: label_bottom_right,
    //};
    //ui.put(label_rect, label);
}


