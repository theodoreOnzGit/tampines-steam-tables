use uom::si::f64::*;

use crate::interfaces::functional_programming::pt_flash_eqm;

/// this is the bread and butter for tampines steam tables, 
/// the control volume
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TampinesSteamTableCV {
    /// these are intensive properties
    pressure: Pressure,
    temperature: ThermodynamicTemperature,
    specific_volume: SpecificVolume,
    specific_enthalpy: AvailableEnergy,
    specific_entropy: SpecificHeatCapacity,

    /// these are extensive properties, 
    /// Volume should be fixed for a control volume by definition
    volume: Volume,
}

impl TampinesSteamTableCV {

    /// creates a new control volume assuming quality is 1 
    /// at the steam table
    ///
    /// this quality is only used at the saturation line of course
    pub fn new_from_tp_quality_1(
        temperature: ThermodynamicTemperature,
        pressure: Pressure,
        volume: Volume,) -> Self {

        let x = 1.0_f64;
        let specific_volume = pt_flash_eqm::v_tp_eqm_two_phase(
            temperature, pressure, x
        );

        let specific_enthalpy = pt_flash_eqm::h_tp_eqm_two_phase(
            temperature, pressure, x
        );

        let specific_entropy = pt_flash_eqm::s_tp_eqm_two_phase(
            temperature, pressure, x
        );

        return Self {
            pressure,
            temperature,
            specific_volume,
            specific_enthalpy,
            specific_entropy,
            volume,
        };
    }
    /// creates a new control volume assuming quality is 0
    /// at the steam table
    ///
    /// this quality is only used at the saturation line of course
    pub fn new_from_tp_quality_0(
        temperature: ThermodynamicTemperature,
        pressure: Pressure,
        volume: Volume,) -> Self {

        let x = 0.0_f64;
        let specific_volume = pt_flash_eqm::v_tp_eqm_two_phase(
            temperature, pressure, x
        );

        let specific_enthalpy = pt_flash_eqm::h_tp_eqm_two_phase(
            temperature, pressure, x
        );

        let specific_entropy = pt_flash_eqm::s_tp_eqm_two_phase(
            temperature, pressure, x
        );

        return Self {
            pressure,
            temperature,
            specific_volume,
            specific_enthalpy,
            specific_entropy,
            volume,
        };
    }



}




