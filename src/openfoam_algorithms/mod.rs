mod simplefoam;

// apparently for two phase flow, 
// this is a useful algorithm to learn from.
//
// https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2008/PraveenPrabhuBaila/Report_twoPhaseEuler.pdf
//
// though for a first model, it may be too complex for a first model
#[allow(non_snake_case)]
mod chtMultiRegionTwoPhaseEulerFoam;


// a less complicated model is drift flux 
// drift flux model
#[allow(non_snake_case)]
mod driftFluxFoam;

// also notable is the VOF methods
// https://www.cfd-online.com/Forums/openfoam-solving/58063-vof-method.html
// which interfoam uses
//
// this helps it to find an interface between the fluid and gas
//
//
// This is interphasechangefoam
// https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2011/MartinAndersen/Tutorial_interPhaseChangeFoam.pdf
//
// apparently accounts for boiling
//
//
// Some other useful tutorials here:
// https://www.wolfdynamics.com/training/mphase/OF2021/mphase_2021_OF8_guided_tutorials.pdf
