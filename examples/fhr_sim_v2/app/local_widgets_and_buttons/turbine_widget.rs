use egui::Widget;

#[derive(Debug,Clone, PartialEq)]
pub struct TurbineWidget {

}

// for the turbine widget, this is basically a collection of rectangles 
// with strokes in there. 
//
// These strokes will vary position with sine wave functions
impl TurbineWidget {

}

impl Widget for TurbineWidget {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {
        todo!()
    }
}
