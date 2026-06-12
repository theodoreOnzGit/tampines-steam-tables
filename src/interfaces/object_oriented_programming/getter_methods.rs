use uom::si::f64::*;
use uom::si::volume::cubic_meter;

use crate::constants::p_crit_water;
use crate::constants::t_crit_water;
use crate::prelude::functional_programming::ph_flash_eqm::ph_flash_region;
use crate::prelude::functional_programming::ph_flash_eqm::x_ph_flash;
use crate::prelude::functional_programming::ph_flash_eqm::w_ph_wood_wallis;
use crate::prelude::functional_programming::ph_flash_eqm::lambda_ph_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::cv_ph_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::cp_ph_eqm;
use crate::dynamic_viscosity::mu_ph_eqm;
use crate::prelude::functional_programming::pt_flash_eqm::FwdEqnRegion;
use crate::region_2_vapour::*;
use crate::region_4_vap_liq_equilibrium::sat_pressure_4;
use crate::region_4_vap_liq_equilibrium::sat_temp_4;
use crate::steam_turbine_equations::choked_flow::single_phase_basic_choked_flow::get_critical_pressure_pure_vapour_ph_stagnation_properties;
use crate::steam_turbine_equations::choked_flow::single_phase_basic_choked_flow::get_critical_pressure_ratio_ideal_gas_using_throat_ph;
use crate::steam_turbine_equations::choked_flow::*;
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

        return w_ph_wood_wallis(p, h);
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
    ///
    /// Note that this is done using throat properties
    pub fn get_critical_pressure_ratio_ideal_gas(&self) -> Ratio {
        
        let p = self.pressure;
        let h = self.specific_enthalpy;

        get_critical_pressure_ratio_ideal_gas_using_throat_ph(p, h)

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

    /// This algorithm uses a more generic approach to 
    /// critical pressure and mass flux,
    ///
    /// basically, one doesn't even find the speed of sound 
    /// but uses a scanning algorithm in order to obtain the 
    /// critical mass flux
    /// this assumes the properties supplied are all stagnation properties
    #[inline]
    pub fn get_crit_pressure_and_massflux(&self) -> (Pressure, MassFlux) {

        let s0 = self.specific_entropy;
        let h0 = self.specific_enthalpy;
        let p0 = self.pressure;

        get_critical_pressure_and_mass_flux_with_stagnation_props(s0, h0, p0)
    }

    /// finds pressure where mach number = 1 during isentropic expansion 
    /// for vapour liquid eqm and subcooled liquid 
    /// it should work vapour as well, just that the vapour algorithm 
    /// tends to use ideal gas critical pressure to bound the search
    /// this one does not
    ///
    ///
    /// Note: I tried mechanical equilibrium where slip ratio = 1 
    /// that is both liquid and vapour move at same velocity 
    /// this does NOT work.
    /// because no matter how low i go in terms of pressure, the 
    /// vle velocity never reaches close to mach 1
    ///
    pub fn get_critical_pressure_vle(&self) -> Pressure {
        let (pressure, _mass_flux) = self.get_crit_pressure_and_massflux();

        return pressure;
    }

    /// Finds the pressure where Mach number = 1 during isentropic expansion
    /// This only works for superheated vapour
    pub fn get_critical_pressure_pure_vapour(&self) -> Pressure {

        let p0 = self.pressure;
        let s0 = self.specific_entropy;
        let h0 = self.specific_enthalpy;

        get_critical_pressure_pure_vapour_ph_stagnation_properties(
            p0, h0, Some(s0)
        )

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
        let h0 = self.get_specific_enthalpy();
        let p0 = self.get_pressure();


        let (_critical_pressure,mass_flux) = 
            get_critical_pressure_and_mass_flux_with_stagnation_props(
                s0, h0, p0
            );


        return mass_flux;
    }
}


