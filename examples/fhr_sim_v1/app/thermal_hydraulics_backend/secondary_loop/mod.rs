use crate::{FHRSimulatorApp};
use crate::app::thermal_hydraulics_backend::fhr_thermal_hydraulics_state::FHRThermalHydraulicsState;
use tampines_steam_tables::interfaces::functional_programming::{ps_flash_eqm, pt_flash_eqm};
use uom::si::f64::*;
use uom::si::pressure::bar;
use uom::si::thermodynamic_temperature::degree_celsius;

impl FHRSimulatorApp {
    /// contains info for calculating secondary loop for one 
    /// timestep 
    ///
    /// the key thing is the tube side temperature
    ///
    /// note that in this simplified steam generator calculation,
    /// everything instantly goes to steady state
    pub(crate) fn secondary_loop_single_timestep(
        fhr_th_state: &FHRThermalHydraulicsState,
        // steam generator settings 
        steam_generator_tube_side_temperature: ThermodynamicTemperature,
        steam_generator_overall_ua: ThermalConductance,
        timestep: Time,
        secondary_loop_mass_flowrate: MassRate,
        // pump settings 
        user_specified_pump_outlet_pressure: Pressure,
    ) -> () {
        let heat_rate_to_steam_generator_tube = 
            -fhr_th_state.heat_added_to_steam_generator_shell_side
            /timestep;
        // first start at the boundary condition
        
        // this is condenser outlet
        // we start at 0.1 bar, 35 degrees C
        let condenser_outlet_pressure = Pressure::new::<bar>(0.1);
        let condenser_outlet_temperature = 
            ThermodynamicTemperature::new::<degree_celsius>(35.0);

        // then calculate the entropy... 
        // 
        let condenser_outlet_entropy: SpecificHeatCapacity
            = pt_flash_eqm::s_tp_eqm_single_phase(
                condenser_outlet_temperature, 
                condenser_outlet_pressure);

        // then condenser goes into pump... with new pressure
        // assume pump is isentropic... for simplicity
        //
        //
        // it pumps it up to minimum of 10 bar 

        let pump_outlet_temperature: ThermodynamicTemperature;

        let ten_bar = Pressure::new::<bar>(10.0);
        let sg_inlet_pressure: Pressure;

        if user_specified_pump_outlet_pressure < ten_bar {
            pump_outlet_temperature = 
                ps_flash_eqm::t_ps_eqm(
                    ten_bar, condenser_outlet_entropy
                );
            sg_inlet_pressure = ten_bar;
        } else {
            pump_outlet_temperature = 
                ps_flash_eqm::t_ps_eqm(
                    user_specified_pump_outlet_pressure, 
                    condenser_outlet_entropy
                );
            sg_inlet_pressure = user_specified_pump_outlet_pressure;

        }
        // now the pump outlet temperature will form the 
        // inlet of the steam generator 

        let steam_generator_tube_inlet_enthalpy: AvailableEnergy 
            = pt_flash_eqm::h_tp_eqm_single_phase(
                pump_outlet_temperature, sg_inlet_pressure
            );

        // for energy balance:
        //
        // mass flowrate * (h_out - h_in) = heat_added_to_steam_gen_tube





    }
}
