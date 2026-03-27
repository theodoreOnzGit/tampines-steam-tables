use uom::si::{f64::*, ratio::ratio};

use crate::{dynamic_viscosity::mu_ph_eqm, prelude::functional_programming::ph_flash_eqm::{cp_ph_eqm, cv_ph_eqm, kappa_ph_eqm, lambda_ph_eqm, w_ph_eqm}};
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

    /// get mach number 
    pub fn get_mach_number(&self, v: Velocity) -> Ratio {

        v/self.get_speed_of_sound()
    }

    /// returns the specific heat ratio cp/cv of steam 
    pub fn get_specific_heat_ratio(&self) -> Ratio {

        let p = self.pressure;
        let h = self.specific_enthalpy;

        let cp = cp_ph_eqm(p, h);
        let cv = cv_ph_eqm(p, h);

        cp/cv
    }

    /// returns cp 
    pub fn get_cp(&self) -> SpecificHeatCapacity {

        let p = self.pressure;
        let h = self.specific_enthalpy;

        cp_ph_eqm(p, h)
    }
    /// returns cv
    pub fn get_cv(&self) -> SpecificHeatCapacity {

        let p = self.pressure;
        let h = self.specific_enthalpy;

        cv_ph_eqm(p, h)
    }

    /// returns thermal thermal_conductivity of steam 
    pub fn get_thermal_conductivity(&self) -> ThermalConductivity {

        let p = self.pressure;
        let h = self.specific_enthalpy;

        lambda_ph_eqm(p, h)
    }

    /// returns critical pressure ratio for choked flow 
    /// ie to accelerate the flow to Mach 1
    ///
    /// P*/P0 = (2/(k+1))^(k/(k-1))
    ///
    /// This is under ideal gas assumption (may not work)
    pub fn get_critical_pressure_ratio_ideal_gas(&self) -> Ratio {
        
        let k = self.get_specific_heat_ratio();

        let ratio_one = Ratio::new::<ratio>(1.0);

        let k_plus_one = k + ratio_one;

        let k_minus_one = k - ratio_one;

        let exponent: f64 = (k/k_minus_one).get::<ratio>();
        let coeff: f64 = (2.0/k_plus_one).get::<ratio>();

        let ratio_value = coeff.powf(exponent);




        Ratio::new::<ratio>(ratio_value)

    }
}

