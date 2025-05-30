// Region 3

const REGION_3_COEFFS: [[f64; 3]; 40] = [
    [0.0, 0.0, 0.10658070028513e1],
    [0.0, 0.0, -0.15732845290239e2],
    [0.0, 1.0, 0.20944396974307e2],
    [0.0, 2.0, -0.76867707878716e1],
    [0.0, 7.0, 0.26185947787954e1],
    [0.0, 10.0, -0.28080781148620e1],
    [0.0, 12.0, 0.12053369696517e1],
    [0.0, 23.0, -0.84566812812502e-2],
    [1.0, 2.0, -0.12654315477714e1],
    [1.0, 6.0, -0.11524407806681e1],
    [1.0, 15.0, 0.88521043984318],
    [1.0, 17.0, -0.64207765181607],
    [2.0, 0.0, 0.38493460186671],
    [2.0, 2.0, -0.85214708824206],
    [2.0, 6.0, 0.48972281541877e1],
    [2.0, 7.0, -0.30502617256965e1],
    [2.0, 22.0, 0.39420536879154e-1],
    [2.0, 26.0, 0.12558408424308],
    [3.0, 0.0, -0.27999329698710],
    [3.0, 2.0, 0.13899799569460e1],
    [3.0, 4.0, -0.20189915023570e1],
    [3.0, 16.0, -0.82147637173963e-2],
    [3.0, 26.0, -0.47596035734923],
    [4.0, 0.0, 0.43984074473500e-1],
    [4.0, 2.0, -0.44476435428739],
    [4.0, 4.0, 0.90572070719733],
    [4.0, 26.0, 0.70522450087967],
    [5.0, 1.0, 0.10770512626332],
    [5.0, 3.0, -0.32913623258954],
    [5.0, 26.0, -0.50871062041158],
    [6.0, 0.0, -0.22175400873096e-1],
    [6.0, 2.0, 0.94260751665092e-1],
    [6.0, 26.0, 0.16436278447961],
    [7.0, 2.0, -0.13503372241348e-1],
    [8.0, 26.0, -0.14834345352472e-1],
    [9.0, 2.0, 0.57922953628084e-3],
    [9.0, 26.0, 0.32308904703711e-2],
    [10.0, 0.0, 0.80964802996215e-4],
    [10.0, 1.0, -0.16557679795037e-3],
    [11.0, 26.0, -0.44923899061815e-4],
];

/// dimensionless temperature (tau)
/// and dimensionless density (delta)
pub mod dimensionless_tau_and_delta;
pub use dimensionless_tau_and_delta::*;

pub mod phi_dimensionless_helmholtz_free_energy;
pub use phi_dimensionless_helmholtz_free_energy::*;


pub mod phi_deriviatives;
pub use phi_deriviatives::*;

/// intensive properties for forward equations in region 3
pub mod intensive_properties;
pub use intensive_properties::*;

/// region_2_3_auxiliary_boundary 
pub mod aux_eqn_boundary_region_2_and_region_3;
pub use aux_eqn_boundary_region_2_and_region_3::*;

/// region 3 ph equations 
pub mod backward_eqn_ph_3;
pub use backward_eqn_ph_3::*;


/// region 3 pt equations for volume 
/// this enables pt flashing in this region 
pub mod backward_eqn_pt_3;
pub use backward_eqn_pt_3::*;

/// region 3 ps equations for volume and temperature
/// this enables pt flashing in this region 
pub mod backward_eqn_ps_3;
pub use backward_eqn_ps_3::*;

/// region 3 hs equations for volume and temperature
/// this enables pt flashing in this region 
pub mod backward_eqn_hs_3;
pub use backward_eqn_hs_3::*;

/// tests 
#[cfg(test)]
mod tests;


