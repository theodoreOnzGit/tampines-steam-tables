//! Test data from Figure 2 of Zaloudek (as reviewed in):
//! Saha, P. (1978). A review of two-phase steam-water critical flow models
//! with emphasis on thermal nonequilibrium. NUREG/CR-0417, BNL-NUREG-50907.
//! Brookhaven National Laboratory, Upton, New York.
//! https://www.nrc.gov/docs/ML1925/ML19256F779.pdf
//!
//! Data format: (critical_pressure_psia, critical_mass_flux_lb_per_s_per_ft2,
//!               stagnation_enthalpy_btu_per_lb)
//! Critical pressures: 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500,
//!                     750, 1000, 1500, 2000, 3000 psia

/// for this test suite, it uses Zaloudek test data, 
/// given the throat pressure and enthalpy, it should guess the stagnation 
/// enthalpy correctly
pub mod backward_throat_to_stagnation;
/// this contains a test suite for stagnation (p0,h0) points within 
/// the vapour liquid equilibrium (VLE)
///
/// this VLE looks like a dome on the (p,h) diagram for the steam 
/// tables.
///
/// Hence this test suite is called in dome stagnation
pub mod in_dome_stagnation;
pub mod generic_multiphase_stagnation;

/// exploratory test suite for subcooled-liquid / liquid-like stagnation
/// points (OUTSIDE the dome, left side), validating the subcooled solver
pub mod subcooled_outside_dome_stagnation;
