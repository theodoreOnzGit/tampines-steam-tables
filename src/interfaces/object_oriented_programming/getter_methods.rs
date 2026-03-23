use uom::si::f64::*;

use crate::{dynamic_viscosity::mu_ph_eqm, prelude::functional_programming::ph_flash_eqm::w_ph_eqm};
impl super::TampinesSteamTableCV {
    /// Returns the pressure of the control volume.
    pub fn get_pressure(&self) -> Pressure {
        self.pressure
    }

    /// Returns the thermodynamic temperature of the control volume.
    pub fn get_temperature(&self) -> ThermodynamicTemperature {
        self.temperature
    }

    /// Returns the specific volume of the fluid in the control volume.
    pub fn get_specific_volume(&self) -> SpecificVolume {
        self.specific_volume
    }

    /// Returns the specific enthalpy of the fluid in the control volume.
    pub fn get_specific_enthalpy(&self) -> AvailableEnergy {
        self.specific_enthalpy
    }

    /// Returns the specific entropy of the fluid in the control volume.
    pub fn get_specific_entropy(&self) -> SpecificHeatCapacity {
        self.specific_entropy
    }

    /// Returns the total volume of the control volume.
    pub fn get_volume(&self) -> Volume {
        self.volume
    }
    /// returns the mass within the control volume 
    pub fn get_mass(&self) -> Mass {
        return self.volume/self.specific_volume;
    }


    /// returns viscosity (important for Reynold's number)
    pub fn get_viscosity(&self) -> DynamicViscosity {

        let p = self.pressure;
        let h = self.specific_enthalpy;

        return mu_ph_eqm(p, h);
    }


    /// returns speed of sound 
    /// important for compressible flow in turbine 
    pub fn get_speed_of_sound(&self) -> Velocity {
        let p = self.pressure;
        let h = self.specific_enthalpy;

        return w_ph_eqm(p, h);
    }
}

