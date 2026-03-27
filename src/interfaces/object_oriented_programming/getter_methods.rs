use uom::{ConstZero, si::{f64::*, pressure::pascal, ratio::ratio}};

use crate::{dynamic_viscosity::mu_ph_eqm, prelude::functional_programming::{ph_flash_eqm::{cp_ph_eqm, cv_ph_eqm, kappa_ph_eqm, lambda_ph_eqm, s_ph_eqm, w_ph_eqm}, ps_flash_eqm::h_ps_eqm}};
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

    /// Returns critical pressure ratio for choked flow using isentropic relations
    /// This accounts for real gas behavior of steam
    ///
    /// The critical pressure P* is found where the flow reaches Mach 1 during
    /// isentropic expansion from stagnation conditions (P0, h0)
    pub fn get_critical_pressure_ratio(&self) -> Ratio {
        let p0 = self.pressure;

        // Find critical pressure where Mach = 1
        // This requires iterative solution
        let p_star = self.find_critical_pressure_isentropic();

        p_star / p0
    }

    /// Finds the pressure where Mach number = 1 during isentropic expansion
    fn find_critical_pressure_isentropic(&self) -> Pressure {

        let ideal_gas_critical_pressure_ratio = 
            self.get_critical_pressure_ratio_ideal_gas();

        let p0 = self.pressure;
        let s0 = self.specific_entropy;
        let h0 = self.specific_enthalpy;
        // Initial guess: use ideal gas approximation as starting point
        let p_guess = p0 * ideal_gas_critical_pressure_ratio; 
        // ~(2/(k+1))^(k/(k-1)) for k≈1.3


        // Newton-Raphson or bisection to find where:
        // v = w (velocity equals speed of sound)
        //
        // From energy equation: h0 = h + v²/2
        // At critical point: v = w, so: h0 = h + w²/2

        let tolerance = Pressure::new::<pascal>(1.0); // 1 Pa tolerance
        let max_iterations = 50;

        // Bisection method bounds
        // Set bounds around the ideal gas guess (±30% to be safe)
        // This reduces iterations compared to starting at 0.1*p0 to 1.0*p0
        let mut p_low = p_guess * 0.7;   // 30% below guess
        let mut p_high = p_guess * 1.3;  // 30% above guess

        // Clamp bounds to reasonable range
        if p_low < p0 * 0.1 {
            p_low = p0 * 0.1;
        }
        if p_high > p0 * 0.99 {
            p_high = p0 * 0.99;
        }


        for _ in 0..max_iterations {
            let p_mid = (p_low + p_high) / 2.0;

            // Get properties at this pressure (isentropic)
            let h_mid = h_ps_eqm(p_mid, s0);
            let w_mid = w_ph_eqm(p_mid, h_mid);

            // Calculate velocity from energy equation
            // h0 = h + v²/2  =>  v = sqrt(2*(h0 - h))
            let delta_h = h0 - h_mid;

            if delta_h < AvailableEnergy::ZERO {
                // Pressure too low, expansion exceeded stagnation enthalpy
                p_low = p_mid;
                continue;
            }

            let v_squared = 2.0 * delta_h;
            let v = v_squared.sqrt();

            // Check if Mach = 1 (v = w)
            let mach = v / w_mid;
            let mach_value = mach.get::<ratio>();

            if (mach_value - 1.0).abs() < 0.0001 {
                return p_mid;
            }

            // Adjust bounds
            if mach_value < 1.0 {
                p_high = p_mid; // Need lower pressure (more expansion)
            } else {
                p_low = p_mid;  // Need higher pressure (less expansion)
            }

            // Check convergence
            if (p_high - p_low) < tolerance {
                return (p_low + p_high) / 2.0;
            }
        }

        // Return midpoint if not converged
        (p_low + p_high) / 2.0
    }
}

