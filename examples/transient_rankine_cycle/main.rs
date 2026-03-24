/// this is a transient rankine cycle, without all the gui 
/// complexities
fn main(){

    isnentropic_nozzle_tests();


}
use tampines_steam_tables::steam_turbine_equations::{force_balance_isentropic_nozzle, get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified};
use uom::si::force::newton;
use uom::si::pressure::bar;
use uom::si::mass_rate::kilogram_per_second;
use uom::si::f64::*;
use uom::si::available_energy::joule_per_kilogram;
use uom::si::area::square_meter;

pub fn isnentropic_nozzle_tests(){

    println!("\n=== HP Turbine First Stage Nozzle ===");

    // Inlet conditions (after superheater)
    let p1 = Pressure::new::<bar>(160.0);  // 160 bar
    let t1_celsius = 540.0;  // 540°C superheated steam

    // For 540°C, 160 bar superheated steam: h ≈ 3,450 kJ/kg
    let h1 = AvailableEnergy::new::<joule_per_kilogram>(3_450_000.0);

    // Mass flow for 250 MW turbine
    let mass_flowrate = MassRate::new::<kilogram_per_second>(200.0);

    // Nozzle areas (converging-diverging)
    let a1 = Area::new::<square_meter>(0.15);  // Inlet area
    let a2 = Area::new::<square_meter>(0.20);  // Outlet area (diverging)

    // now let's plot a trend of isentropic nozzles 

    println!("{:?}",&("oulet pressure (Bar)","force bal(Newton)"));
    let mut p2 = p1;
    let n = 100;
    for i in 0..n {

        p2 = (n as f64 - i as f64)/(n as f64) * p1;

        let force_bal: Force = 
            force_balance_isentropic_nozzle(p1, p2, h1, mass_flowrate, a1, a2);

        println!("{:?}",&(p2.get::<bar>(),force_bal.get::<newton>()));

    }

    // let's also check the lower bound pressures
    // below the first 100 points
    for i in 0..n {

        p2 = (1_f64)/(n as f64) * p1;
        p2 *= (n-i) as f64 / (n as f64);


        let force_bal: Force = 
            force_balance_isentropic_nozzle(p1, p2, h1, mass_flowrate, a1, a2);

        println!("{:?}",&(p2.get::<bar>(),force_bal.get::<newton>()));

    }

    //let p2 = get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified(
    //    p1, h1, mass_flowrate, a1, a2, Option::None
    //);
}


// this is vibe coded 

use tampines_steam_tables::prelude::TampinesSteamTableCV;
use uom::si::f64::*;
// ... other imports

struct RankineSystem {
    pump: TampinesSteamTableCV,
    boiler: TampinesSteamTableCV,
    turbine: TampinesSteamTableCV,
    condenser: TampinesSteamTableCV,
    
    // System parameters
    pump_curve: fn(MassRate) -> Pressure, // A function describing the pump's head
    friction_factors: (f64, f64, f64), // For boiler, turbine, condenser
    heat_source_power: Power, // e.g., in Watts
}

struct SolverControls {
    delta_t: Time,
    end_time: Time,
    n_outer_correctors: i32, // PIMPLE loops
    n_inner_correctors: i32, // PISO loops
}
