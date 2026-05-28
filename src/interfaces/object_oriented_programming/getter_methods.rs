use uom::ConstZero;
use uom::si::f64::*;
use uom::si::pressure::pascal;
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::constants::p_crit_water;
use crate::constants::t_crit_water;
use crate::interfaces::functional_programming::hs_flash_eqm::p_hs_eqm;
use crate::interfaces::functional_programming::ps_flash_eqm::w_ps_eqm;
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
    /// This only works for vapour
    pub fn get_critical_pressure_ratio_pure_vapour(&self) -> Ratio {
        let p0 = self.pressure;

        // Find critical pressure where Mach = 1
        // This requires iterative solution
        let p_star = self.get_critical_pressure_pure_vapour();

        p_star / p0
    }

    /// finds pressure where mach number = 1 during isentropic expansion 
    /// for vapour liquid eqm and subcooled liquid 
    /// it should work vapour as well, just that the vapour algorithm 
    /// tends to use ideal gas critical pressure to bound the search
    /// this one does not
    pub fn get_critical_pressure_vle(&self) -> Pressure {

        // for this, the same thing applies 
        // we have isentropic expansion (ie reduction of pressure)
        // such that the mach value is 1
        //

        // first we get stagnation properties 
        // and p0 will be the high bound pressure
        let p0 = self.pressure;
        let s0 = self.specific_entropy;
        let h0 = self.specific_enthalpy;

        let mut p_high = p0;
        let mut p_low = 0.1 * p_high;

        // at stagnation pressure, the pressure would be the highest 
        // so it is closest to liquid 
        // so the speed of sound is the highest 

        let root_finder_pressure = |p_test: Pressure| -> f64 {

            let h_test = h_ps_eqm(p_test, s0);
            let w_test = w_ph_eqm(p_test, h_test);
            // Calculate velocity from energy equation
            // h0 = h + v²/2  =>  v = sqrt(2*(h0 - h))
            let delta_h = h0 - h_test;

            let v_squared = 2.0 * delta_h;
            let v = v_squared.sqrt();
            // Check if Mach = 1 (v = w)
            let mach = v / w_test;
            let mach_value = mach.get::<ratio>();

            return mach_value - 1.0;
        };

        // the high bound for velocity is the speed of sound at stagnation,
        // which should be the highest possible
        // I am giving a 30% factor up
        // the lowest is velocity = 0 m/s (stagnation)
        let mut v_upper_limit = self.get_speed_of_sound() * 1.3;
        let mut v_lower_limit = Velocity::ZERO;
        let v_decrement = v_upper_limit * 0.05;
        let mut v_test = v_upper_limit - v_decrement;

        let root_finder_velocity = |v_test: Velocity| -> f64 {
            let h_test = h0 - 0.5 * v_test * v_test;
            let p_test = p_hs_eqm(h_test, s0);
            let w_test = w_ps_eqm(p_test, s0);


            // Check if Mach = 1 (v = w)
            let mach = v_test / w_test;
            let mach_value = mach.get::<ratio>();

            return mach_value - 1.0;
        };
        // we shall test for sign change 

        let mach_error_initial: f64 = root_finder_velocity(v_test);

        let debug = false;

        // if i were to use a velocity scanner, I would go from supersonic 
        // speed down to subsonic
        //
        // here is my velocity scanner

        // so I'm going to do the velocity scan first

        while v_test > v_lower_limit {

            // we are going down to mach 1

            
            let mach_error: f64 = root_finder_velocity(v_test);

            // check if signs are same, then continue
            if mach_error * mach_error_initial >= 0.0 {
                v_upper_limit = v_test;
                v_test -= v_decrement;

                if debug {
                    dbg!(&(v_test));
                    dbg!(&(mach_error));
                }
                continue;
            }

            // if signs are not the same, then break out

            // if signs are not same

            if debug {
                dbg!(&(v_test));
                dbg!(&(mach_error));
            }
            v_lower_limit = v_test; 
            break;
        };
        
        // now i can do bisection (or a secant method) 
        // between these two limits
        // since it's quite near the root
        //
        // or as AI suggested, I'm going to try Regula Falsi
        // near this region

        if debug{
            println!("Regula Falsi bounds found");
            dbg!(&(v_lower_limit,v_upper_limit));
        }

        let tolerance = Pressure::new::<pascal>(1.0); // 1 Pa tolerance
        let max_iterations = 50;

        let mut error_lower_limit = root_finder_velocity(v_lower_limit);
        let mut error_upper_limit = root_finder_velocity(v_upper_limit);
        // this time i use regula falsi
        //
        // this ensures the bounds are not the same sign 
        // just a sanity check
        if error_lower_limit * error_upper_limit>= 0.0 {
            panic!("bounds are same sign!");
        }

        const TOLERANCE: f64 = 0.0001;  // 0.01% tolerance
        // this is regula falsi
        for _ in 0..max_iterations {

            // using secant formula
            v_test = 
                v_upper_limit - (v_upper_limit - v_lower_limit) * 
                error_upper_limit/(error_upper_limit - error_lower_limit);
            // check mach number error 

            let mach_error = root_finder_velocity(v_test);

            if mach_error.abs() < TOLERANCE {
                // we found the critical pressure 
                let h_test = h0 - 0.5 * v_test * v_test;
                let p_test = p_hs_eqm(h_test, s0);
                return p_test;
            }

            // if not, keep updating the bounds
            // keep root bracketed 

            if error_lower_limit * mach_error < 0.0 {

                v_upper_limit = v_test;
                error_upper_limit = mach_error;
            } else {

                v_lower_limit = v_test;
                error_lower_limit = mach_error;
            }


            if debug {
                dbg!(&(v_lower_limit,v_upper_limit));
            }


        }
        panic!("unable to find critical pressure");

    }

    /// Finds the pressure where Mach number = 1 during isentropic expansion
    /// This only works for vapour
    pub fn get_critical_pressure_pure_vapour(&self) -> Pressure {

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

    /// critical mass flux 
    /// for choked flow
    /// assumes state supplied is stagnation state
    pub fn get_stagnation_critical_mass_flux(&self) -> MassFlux {


        let s0 = self.get_specific_entropy();
        let s1 = s0;


        // now, i'll have to get a solver for choked flow 

        // let's use the critical pressure 


        // this is critical pressure for mach 1
        let p2 = self.get_critical_pressure_pure_vapour();
        // let's get speed of sound here 
        let s2 = s1;
        let v2 = self.get_volume();
        let state_2 = Self::new_from_ps(p2, s2, v2);
        let c = state_2.get_speed_of_sound();
        let rho_2 = state_2.get_rho();

        return c*rho_2;
    }
}


