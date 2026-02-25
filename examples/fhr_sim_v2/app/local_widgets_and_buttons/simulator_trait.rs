use egui::Vec2;
pub trait VisualSimulationObject {


    /// this function is here to ensure that the coordinates at the  
    /// end joint are given
    fn get_end_joint_coordinate(&self) -> Vec2;
}
