use uom::si::{f64::*, mass::kilogram};
impl super::TampinesSteamTableCV {
    pub fn set_tpx(&mut self, t: ThermodynamicTemperature, p: Pressure, x: f64) {
        let volume = self.get_volume();

        *self = Self::new_from_tp_quality(t, p, volume, x);
    }

    pub fn set_ph(&mut self, p: Pressure, h: AvailableEnergy) {
        let volume = self.get_volume();

        *self = Self::new_from_ph(p, h, volume);
    }
    pub fn set_ps(&mut self, p: Pressure, s: SpecificHeatCapacity) {
        let volume = self.get_volume();

        *self = Self::new_from_ps(p, s, volume);
    }

    // AI assistance:
    // Add these to your setter methods implementation block

    /// Models an ideal (isentropic) compression process to a new pressure.
    /// The entropy is held constant. Perfect for an ideal pump.
    pub fn compress_isentropically(&mut self, new_p: Pressure) {
        let current_s = self.get_specific_entropy();
        self.set_ps(new_p, current_s);
    }

    /// Models an ideal (isentropic) expansion process to a new pressure.
    /// The entropy is held constant. Perfect for an ideal turbine.
    pub fn expand_isentropically(&mut self, new_p: Pressure) {
        let current_s = self.get_specific_entropy();
        self.set_ps(new_p, current_s);
    }

    /// Models a process of adding heat at constant pressure (isobaric).
    /// This is ideal for a boiler.
    pub fn add_heat_isobaric(&mut self, heat_added_per_kg: AvailableEnergy) {
        let current_p = self.get_pressure();
        let new_h = self.get_specific_enthalpy() + heat_added_per_kg;
        self.set_ph(current_p, new_h);
    }

    /// Models a process of removing heat at constant pressure (isobaric).
    /// This is ideal for a condenser.
    pub fn remove_heat_isobaric(&mut self, heat_removed_per_kg: AvailableEnergy) {
        let current_p = self.get_pressure();
        let new_h = self.get_specific_enthalpy() - heat_removed_per_kg;
        self.set_ph(current_p, new_h);
    }

    // AI assistance:
    //
    // --- New Extensive Process Methods ---

    /// Models adding a total amount of heat (extensive) to the control volume
    /// at constant pressure (isobaric process).
    ///
    /// # Arguments
    /// * `total_heat_added`: The total energy added to the control volume, in Joules.
    pub fn add_heat_isobaric_extensive(&mut self, total_heat_added: Energy) {
        let mass = self.get_mass();

        // Safety check to prevent division by zero if the control volume is empty
        if mass <= Mass::new::<kilogram>(0.0) {
            return;
        }

        // Convert extensive heat (J) to intensive heat (J/kg)
        let heat_added_per_kg = total_heat_added / mass;

        let current_p = self.get_pressure();
        let new_h = self.get_specific_enthalpy() + heat_added_per_kg;

        // Use the existing setter to update the state
        self.set_ph(current_p, new_h);
    }

    /// Models removing a total amount of heat (extensive) from the control volume
    /// at constant pressure (isobaric process).
    ///
    /// # Arguments
    /// * `total_heat_removed`: The total energy removed from the control volume, in Joules.
    pub fn remove_heat_isobaric_extensive(&mut self, total_heat_removed: Energy) {
        let mass = self.get_mass();

        if mass <= Mass::new::<kilogram>(0.0) {
            return;
        }

        // Convert extensive heat (J) to intensive heat (J/kg)
        let heat_removed_per_kg = total_heat_removed / mass;

        let current_p = self.get_pressure();
        let new_h = self.get_specific_enthalpy() - heat_removed_per_kg;

        self.set_ph(current_p, new_h);
    }

    /// Models the work done on the fluid during compression by specifying the
    /// total work input and the target pressure. This is ideal for a real pump.
    ///
    /// # Arguments
    /// * `total_work_input`: The total work done ON the fluid, in Joules.
    /// * `outlet_pressure`: The target pressure after compression.
    pub fn compress_with_work_extensive(
        &mut self,
        total_work_input: Energy,
        outlet_pressure: Pressure,
    ) {
        let mass = self.get_mass();

        if mass <= Mass::new::<kilogram>(0.0) {
            return;
        }

        // Convert extensive work (J) to intensive work (J/kg)
        let work_input_per_kg = total_work_input / mass;

        // The first law for a control volume: w_in = h_out - h_in
        let new_h = self.get_specific_enthalpy() + work_input_per_kg;

        // The new state is defined by the new enthalpy and the target outlet pressure
        self.set_ph(outlet_pressure, new_h);
    }

    /// Models the work done by the fluid during expansion by specifying the
    /// total work output and the target pressure. This is ideal for a real turbine.
    ///
    /// # Arguments
    /// * `total_work_output`: The total work done BY the fluid, in Joules.
    /// * `outlet_pressure`: The target pressure after expansion.
    pub fn expand_with_work_extensive(
        &mut self,
        total_work_output: Energy,
        outlet_pressure: Pressure,
    ) {
        let mass = self.get_mass();

        if mass <= Mass::new::<kilogram>(0.0) {
            return;
        }

        // Convert extensive work (J) to intensive work (J/kg)
        let work_output_per_kg = total_work_output / mass;

        // The first law for a control volume: w_out = h_in - h_out
        let new_h = self.get_specific_enthalpy() - work_output_per_kg;

        // The new state is defined by the new enthalpy and the target outlet pressure
        self.set_ph(outlet_pressure, new_h);
    }
}
