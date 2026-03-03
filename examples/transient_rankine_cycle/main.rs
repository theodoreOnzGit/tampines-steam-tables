/// this is a transient rankine cycle, without all the gui 
/// complexities
fn main(){

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
