use uom::si::available_energy::joule_per_kilogram;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::{f64::*, ratio::ratio};
use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};

/// this struct helps to manage pebble bed thermal hydraulics 
/// calculations
///
/// This is important for fast feedback of fuel temperature 
///
/// basically it keeps track of the thermal inertia of the system
/// among other things
///
#[derive(Clone, Debug, Copy)]
pub struct PebbleBedThermalHydraulics {
    /// this is the current specific enthalpy of the pebble bed
    pub current_fuel_specific_enthalpy: AvailableEnergy,
}

impl PebbleBedThermalHydraulics {

    /// simple function to calculate enthalpy at next timestep 
    /// and heat removal from pebble bed
    pub fn calc_th_and_return_heat_removal_from_pebble_bed(
        &mut self, 
        timestep: Time,
        fission_power_adjusted_for_decay_heat: Power,
        pebble_bed_mass: Mass,
        pebble_bed_heat_transfer_area: Area,
        pebble_bed_overall_htc: HeatTransfer,
        pebble_bed_coolant_temp: ThermodynamicTemperature,
        ) -> Power {

        let heat_rate_transferred_from_pebble_bed_to_coolant: Power = 
            - ( pebble_bed_overall_htc * pebble_bed_heat_transfer_area * pebble_bed_coolant_temp 
                - pebble_bed_overall_htc * pebble_bed_heat_transfer_area * self.get_temperature_from_enthalpy_uo2_heuristic()
            );


        let net_power_gain = 
            fission_power_adjusted_for_decay_heat 
            - heat_rate_transferred_from_pebble_bed_to_coolant;

        let new_pebble_bed_enthalpy = 
            self.current_fuel_specific_enthalpy * pebble_bed_mass 
            + net_power_gain * timestep;

        let new_pebble_bed_specific_entahlpy = 
            new_pebble_bed_enthalpy/pebble_bed_mass;

        self.current_fuel_specific_enthalpy = new_pebble_bed_specific_entahlpy;

        return heat_rate_transferred_from_pebble_bed_to_coolant;
        
    }

    /// creates a new instance of PebbleBedThermalHydraulics
    pub fn new() -> Self {
        // let's start pebble bed at 500 degc 
        let start_temp = ThermodynamicTemperature::new::<degree_celsius>(500.0);
        let start_enthalpy = Self::get_enthalpy_from_temperature_uo2(start_temp);

        return Self { current_fuel_specific_enthalpy: start_enthalpy };
    }

    /// based on assumed constant cp of 340 J/(kg K)
    /// enthalpy at 298K is deemd to be zero J/kg
    pub fn get_enthalpy_from_temperature_uo2_heuristic(
        fuel_temp: ThermodynamicTemperature,
    ) -> AvailableEnergy {

        let ref_temp_298 = TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(298.0);
        let temp_minus_298 = fuel_temp - ref_temp_298;
        let cp = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(340.0);
        return temp_minus_298 * cp;

    }
    /// based on assumed constant cp of 340 J/(kg K)
    /// enthalpy at 298K is deemd to be zero J/kg
    pub fn get_temperature_from_enthalpy_uo2_heuristic(
        &self
    ) -> ThermodynamicTemperature {

        // h = cp (T - 298) 
        let ref_temp_298 = TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(298.0);
        let cp = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(340.0);
        let cp_times298: AvailableEnergy = cp * ref_temp_298;

        let h_plus_cp_times_298: AvailableEnergy = self.current_fuel_specific_enthalpy
            + cp_times298;

        let fuel_temp_estimate: ThermodynamicTemperature = 
            ThermodynamicTemperature::new::<kelvin>(
                (h_plus_cp_times_298/cp).get::<uom::si::temperature_interval::kelvin>()
            );

        return fuel_temp_estimate;


    }

    /// Carbajo, J. J. (2001). Thermophysical properties of MOX and UO2 
    /// fuels including the effects of irradiation (No. ORNL/TM-2000/351). 
    /// Oak Ridge National Lab.(ORNL), Oak Ridge, TN (United States).
    ///
    /// enthalpy at 298K is deemd to be zero J/kg
    pub fn get_enthalpy_from_temperature_uo2(
        fuel_temp: ThermodynamicTemperature
    ) -> AvailableEnergy {

        // see page 18
        // table 4.2
        // for uranium oxide

        let c1 = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
            302.27
        );
        let c2 = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
            8.463e-3
        ) / TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(1.0);
        let c3 = AvailableEnergy::new::<joule_per_kilogram>(
            8.741e7
        );

        // einstein temp
        let theta = ThermodynamicTemperature::new::<kelvin>(
            548.68
        );

        // electron activation energy divide by boltzmann constant 
        let e_a = ThermodynamicTemperature::new::<kelvin>(
            18_531.7
        );

        // some terms for calc 

        let ref_temp_298 = ThermodynamicTemperature::new::<kelvin>(298.0);
        let theta_by_t: Ratio = theta/fuel_temp;
        let theta_by_298: Ratio = theta/ref_temp_298;

        let e_a_by_t: Ratio = e_a/fuel_temp;

        // calc term by term
        let term1: AvailableEnergy = c1 * theta * (
            (theta_by_t.get::<ratio>().exp() - 1.0).recip()
            -(theta_by_298.get::<ratio>().exp() - 1.0).recip()
        );

        let term2: AvailableEnergy = 
            c2 * (fuel_temp * fuel_temp - ref_temp_298 * ref_temp_298);

        let term3: AvailableEnergy = c3 * (-e_a_by_t).exp();

        // enthalpy increment over 298K 
        let enthalpy_increment_over_298_kelvin = 
            term1 + term2 + term3;

        return enthalpy_increment_over_298_kelvin;




    }

    /// assumes uo2 in pebble, gets temp
    /// based on some backward relation (no iteration)
    pub fn get_current_temeperature_uo2(&self) -> ThermodynamicTemperature {
        todo!()
    }
}

