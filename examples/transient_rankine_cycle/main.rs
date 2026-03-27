/// this is a transient rankine cycle, without all the gui 
/// complexities
fn main(){

    isnentropic_nozzle_tests();


}
use tampines_steam_tables::prelude::functional_programming::ph_flash_eqm::s_ph_eqm;
use tampines_steam_tables::prelude::functional_programming::ps_flash_eqm::h_ps_eqm;
use tampines_steam_tables::steam_turbine_equations::{force_balance_isentropic_nozzle, get_isentropic_nozzles_outlet_ph_rho_point_ps_algo_simplified};
use uom::si::force::newton;
use uom::si::pressure::{bar, pascal};
use uom::si::mass_rate::kilogram_per_second;
use uom::si::f64::*;
use uom::si::available_energy::joule_per_kilogram;
use uom::si::area::square_meter;

pub fn isnentropic_nozzle_tests(){

}

