/// This is a 3 phase generator where 
/// phase shift is 0 degrees, 60 degrees and 120 degrees 
/// respectively
#[allow(non_snake_case)]
use uom::si::f64::*;
#[derive(Debug,Clone, PartialEq)]
pub struct ThreePhaseElectricGenerator {
    /// I
    I: MomentOfInertia,

    /// omega 
    omega: AngularVelocity,

    /// no of turns for the coil
    N: usize,

    /// magnetic flux density for generator
    B: MagneticFluxDensity,

    /// area of coil 
    A: Area,

    /// turbine efficiency 
    eta: Ratio,


}

impl ThreePhaseElectricGenerator {

    /// this immutably calculates new angular velocity 
    /// in an explicit manner, given a source term
    pub fn calculate_new_angular_velocity(
        &self,
        source: Torque,
        load_resistance: ElectricalResistance,
        current_time: Time,
        delta_t: Time,
    ) -> AngularVelocity {

        todo!();
    }

}
