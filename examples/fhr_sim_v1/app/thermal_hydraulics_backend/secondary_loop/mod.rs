use crate::{FHRSimulatorApp};
use crate::app::thermal_hydraulics_backend::fhr_thermal_hydraulics_state::FHRThermalHydraulicsState;
use tampines_steam_tables::interfaces::functional_programming::{ph_flash_eqm, ps_flash_eqm, pt_flash_eqm};
use uom::si::f64::*;
use uom::si::mass_rate::kilogram_per_second;
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
        fhr_th_state: &mut FHRThermalHydraulicsState,
        timestep: Time,
        user_specified_secondary_loop_mass_flowrate: &mut MassRate,
        // pump settings 
        user_specified_pump_outlet_pressure: Pressure,
    ) -> SecondaryLoopState {
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
        //
        // now, if mass flowrate is less than some value, crash the 
        // system because the steam generator is too hot 
        if *user_specified_secondary_loop_mass_flowrate <
            MassRate::new::<kilogram_per_second>(0.1) {
                println!("steam generator feedwater flowrate too low!");
                panic!("mass flowrate <0.1 kg/s, unsafe operation");
        }
        //
        let steam_generator_tube_outlet_enthalpy = 
            heat_rate_to_steam_generator_tube/
            *user_specified_secondary_loop_mass_flowrate 
            + steam_generator_tube_inlet_enthalpy;
        if *user_specified_secondary_loop_mass_flowrate <
            MassRate::new::<kilogram_per_second>(0.1) {
                println!("steam generator feedwater flowrate too low!");
                panic!("mass flowrate <0.1 kg/s, unsafe operation");
        }

        let sg_outlet_pressure = sg_inlet_pressure;

        let steam_gen_tube_outlet_temperature = 
            ph_flash_eqm::t_ph_eqm(
                sg_outlet_pressure, steam_generator_tube_outlet_enthalpy
            );

        // then with this enthalpy and pressure, have an isentropic 
        // enthalpy drop in the turbine 
        // to specified pressure (0.2 bar) 

        let turbine_inlet_enthalpy = steam_generator_tube_outlet_enthalpy;
        let turbine_inlet_pressure = sg_outlet_pressure;
        let turbine_inlet_entropy = ph_flash_eqm::s_ph_eqm(
            turbine_inlet_pressure, steam_generator_tube_outlet_enthalpy
        );

        // this is fixed... can modify in future to make this user changeable
        let turbine_outlet_pressure = 
            Pressure::new::<bar>(0.5);
        let turbine_outlet_entropy = turbine_inlet_entropy;
        

        let turbine_outlet_enthalpy = 
            ps_flash_eqm::h_ps_eqm(
                turbine_outlet_pressure, 
                turbine_outlet_entropy
            );

        let work_done_by_turbine_per_unit_mass = 
            turbine_inlet_enthalpy - turbine_outlet_enthalpy;
        let turbine_power: Power = 
            *user_specified_secondary_loop_mass_flowrate * work_done_by_turbine_per_unit_mass;

        // now it goes into the condenser 
        let condenser_inlet_enthalpy = turbine_outlet_enthalpy;
        let condenser_outlet_enthalpy = 
            pt_flash_eqm::h_tp_eqm_single_phase(
                condenser_outlet_temperature, condenser_outlet_pressure
            );

        // this is the required cooling rate for condenser at steady state
        let condenser_duty: Power = 
            *user_specified_secondary_loop_mass_flowrate * 
            (condenser_inlet_enthalpy - condenser_outlet_enthalpy);


        let secondary_loop_state = 
            SecondaryLoopState {
                steam_gen_tube_outlet_temperature,
                turbine_power,
                condenser_duty,
            };


        return secondary_loop_state;

    }
}

/// contains important info for various temperatures and pressures 
/// in the secondary loop
#[derive(Debug,Clone)]
pub struct SecondaryLoopState {
    /// steam generator tube outlet temperature (boiler temperature)
    pub steam_gen_tube_outlet_temperature: ThermodynamicTemperature,
    /// turbine power output 
    pub turbine_power: Power,
    /// condenser cooling duty at steady state
    pub condenser_duty: Power,
}
