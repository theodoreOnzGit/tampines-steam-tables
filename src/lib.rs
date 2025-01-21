pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}

/// constants for the steam table calculations
pub mod constants;

/// region 1 
///
/// Temperature from 273.15 to 623.15 K 
/// Pressure from 0 to 100 MPa
///
/// Up to the saturation line. 
/// This I believe is subcooled liquid region
pub mod region_1_subcooled_liquid;



/// region 2 
///
/// Saturated water region
pub mod region_2_saturated_water;

/// region 3 
///
/// vapour liquid equilibrium
/// wet steam region 
///
/// auxilliary equation for region 2 and 3 are also put here
pub mod region_3_vap_liq_mixture;

/// region 4
///
/// dry saturated steam region 
pub mod region_4_dry_saturated_steam;


/// region 5 
///
/// superheated steam region 
pub mod region_5_superheated_steam;
