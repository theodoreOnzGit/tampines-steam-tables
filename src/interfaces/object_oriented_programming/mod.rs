use uom::si::f64::*;

use crate::{interfaces::functional_programming::*, prelude::functional_programming::hs_flash_eqm::p_hs_eqm, region_4_vap_liq_equilibrium::{sat_pressure_4, sat_temp_4}};

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
    pub fn new_from_tp_quality(
        temperature: ThermodynamicTemperature,
        pressure: Pressure,
        volume: Volume,
        x: f64) -> Self {

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


    pub fn new_from_ph(
        p: Pressure,
        h: AvailableEnergy,
        volume: Volume) -> Self {

        let t = ph_flash_eqm::t_ph_eqm(p, h);
        let specific_volume = ph_flash_eqm::v_ph_eqm(p, h);
        let pressure = p;
        let temperature = t;
        let specific_enthalpy = h;
        let specific_entropy = ph_flash_eqm::s_ph_eqm(p, h);

        return Self {
            pressure,
            temperature,
            specific_volume,
            specific_enthalpy,
            specific_entropy,
            volume,
        };

    }

    pub fn new_from_ps(
        p: Pressure,
        s: SpecificHeatCapacity,
        volume: Volume) -> Self {

        let t = ps_flash_eqm::t_ps_eqm(p, s);
        let specific_volume = ps_flash_eqm::v_ps_eqm(p, s);
        let pressure = p;
        let temperature = t;
        let specific_enthalpy = ps_flash_eqm::h_ps_eqm(p, s);
        let specific_entropy = s;

        return Self {
            pressure,
            temperature,
            specific_volume,
            specific_enthalpy,
            specific_entropy,
            volume,
        };

    }


    pub fn new_from_hs(
        h: AvailableEnergy,
        s: SpecificHeatCapacity,
        volume: Volume) -> Self {

        let p = p_hs_eqm(h, s);
        return Self::new_from_ph(p, h, volume);
    }

    pub fn new_from_sat_pressure_quality(
        p: Pressure,
        x: f64,
        volume: Volume,) -> Self {

        let t: ThermodynamicTemperature = sat_temp_4(p);

        return Self::new_from_tp_quality(
            t, p, volume, x
        );

    }

}




/// vibe coded getter methods
pub mod getter_methods;

/// setter methods 
/// this will deal with setting new thermodynamic equilibrium
/// based on user input parameters
pub mod setter_methods;
