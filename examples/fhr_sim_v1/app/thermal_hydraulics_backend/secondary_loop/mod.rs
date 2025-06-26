use crate::{FHRSimulatorApp};
use crate::app::thermal_hydraulics_backend::fhr_thermal_hydraulics_state::FHRThermalHydraulicsState;
use uom::si::f64::*;

impl FHRSimulatorApp {
    /// contains info for calculating secondary loop for one 
    /// timestep 
    ///
    /// the key thing is the tube side temperature
    pub(crate) fn secondary_loop_single_timestep(
        fhr_th_state: &FHRThermalHydraulicsState,
        // steam generator settings 
        steam_generator_tube_side_temperature: ThermodynamicTemperature,
        steam_generator_overall_ua: ThermalConductance,
    ) -> () {
        // first start at the boundary condition
    }
}
