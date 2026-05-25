use uom::ConstZero;
use uom::si::f64::*;
use uom::si::pressure::pascal;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::constants::p_crit_water;
use crate::constants::t_crit_water;
use crate::prelude::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::prelude::functional_programming::ph_flash_eqm::x_ph_flash;
use crate::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::w_ph_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::lambda_ph_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::cv_ph_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::cp_ph_eqm;
use crate::dynamic_viscosity::mu_ph_eqm;
use crate::prelude::functional_programming::pt_flash_eqm::FwdEqnRegion;
use crate::region_2_vapour::*;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
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
        let p_star = self.get_critical_pressure();

        p_star / p0
    }

    /// Finds the pressure where Mach number = 1 during isentropic expansion
    pub fn get_critical_pressure(&self) -> Pressure {

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

    pub fn get_rho(&self) -> MassDensity {
        self.get_specific_volume().recip()
    }

    // get region of steam 
    pub fn get_region(&self) -> FwdEqnRegion {
        
        let p = self.pressure;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        return region;

    }

    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_specific_volume(&self) -> Option<SpecificVolume>{

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let v = v_tp_2_metastable(t, p);
                return Some(v);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }

    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_specific_enthalpy(&self) -> 
        Option<AvailableEnergy>
    {

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let h = h_tp_2_metastable(t, p);
                return Some(h);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }


    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_internal_energy(&self) -> 
        Option<AvailableEnergy>
    {

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let u = u_tp_2_metastable(t, p);
                return Some(u);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }


    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_specific_entropy(&self) -> 
        Option<SpecificHeatCapacity>
    {

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let s = s_tp_2_metastable(t, p);
                return Some(s);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }
    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_cp(&self) -> 
        Option<SpecificHeatCapacity>
    {

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let cp = cp_tp_2_metastable(t, p);
                return Some(cp);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }
    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_cv(&self) -> 
        Option<SpecificHeatCapacity>
    {

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let cv = cv_tp_2_metastable(t, p);
                return Some(cv);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }
    /// get metastable steam state, (region 2 only) 
    /// 
    /// if not region 2, then returns a None value
    pub fn get_metastable_steam_speed_of_sound(&self) -> 
        Option<Velocity>
    {

        let p = self.pressure;
        let t = self.temperature;
        let h = self.specific_enthalpy;
        let region = ph_flash_region(p, h);

        match region {
            FwdEqnRegion::Region2 => {
                let c = w_tp_2_metastable(t, p);
                return Some(c);
            },
            FwdEqnRegion::Region1 => None,
            FwdEqnRegion::Region3 => None,
            FwdEqnRegion::Region4 => None,
            FwdEqnRegion::Region5 => None,
        }
    }

    
    /// get the steam quality, only if the region is in region 4
    /// region 4 is the vapour liquid equilibrium
    pub fn get_quality(&self) -> f64{

        let p = self.pressure;
        let h = self.specific_enthalpy;

        let x = x_ph_flash(p,h);
        x
    }
    /// get the saturation temperature based on pressure 
    /// provided pressure is less than p_crit
    pub fn try_new_tsat_based_on_pressure(&self) -> Option<ThermodynamicTemperature>{
        let p_crit = p_crit_water();

        if self.pressure > p_crit {
            return None;
        }

        if self.pressure == p_crit {
            return Some(t_crit_water());
        }

        let tsat = sat_temp_4(self.pressure);

        return Some(tsat);
    }

    /// get the saturation pressure based on temperature 
    /// provided temperature is less than t_crit
    pub fn try_new_psat_based_on_temperature(&self) -> Option<Pressure>{
        let t_crit = t_crit_water();

        if self.temperature > t_crit {
            return None;
        }

        if self.temperature == t_crit {
            return Some(p_crit_water());
        }

        let psat = sat_pressure_4(self.temperature);

        return Some(psat);
    }

    /// get the saturation temperature based on pressure 
    /// provided pressure is less than p_crit
    pub fn try_get_tsat(p: Pressure) -> Option<ThermodynamicTemperature>{
        let p_crit = p_crit_water();

        if p > p_crit {
            return None;
        }

        if p == p_crit {
            return Some(t_crit_water());
        }

        let tsat = sat_temp_4(p);

        return Some(tsat);
    }

    /// get the saturation pressure based on temperature 
    /// provided temperature is less than t_crit
    pub fn try_get_psat(t: ThermodynamicTemperature) -> Option<Pressure>{
        let t_crit = t_crit_water();

        if t > t_crit {
            return None;
        }

        if t == t_crit {
            return Some(p_crit_water());
        }

        let psat = sat_pressure_4(t);

        return Some(psat);
    }

    /// just a convenience function to get ref volume 
    /// 1m3
    pub fn get_ref_vol() -> Volume {
        Volume::new::<cubic_meter>(1.0)
    }
}


