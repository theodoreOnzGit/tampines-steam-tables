//// Solve the Momentum equation
//
//MRF.correctBoundaryVelocity(U);
//
//tmp<fvVectorMatrix> tUEqn
//(
//    fvm::ddt(rho, U) + fvm::div(phi, U)
//  + MRF.DDt(rho, U)
//  + turbulence->divDevRhoReff(U)
// ==
//    fvOptions(rho, U)
//);
//fvVectorMatrix& UEqn = tUEqn.ref();
//
//UEqn.relax();
//
//fvOptions.constrain(UEqn);
//
//if (pimple.momentumPredictor())
//{
//    solve(UEqn == -fvc::grad(p));
//
//    fvOptions.correct(U);
//    K = 0.5*magSqr(U);
//}
//
use uom::si::f64::*;

/// within UEqn, we work with mass flowrates
/// note: this wasn't vibe coded
#[derive(Debug, Clone, PartialEq)]
pub struct UEqn {
    /// note that for H, we need 
    /// H = m^n/delta_t - rho g A_{XS} dz/dx
    /// we note that dz is the elevation vector 
    /// dz = sin theta dx
    delta_t: Time,
    mass_flowrate_vector_last_iter: Vec<MassRate>,
    density_vector_last_iter: Vec<MassDensity>,
    g: Acceleration,
    cross_sectional_area_vec_last_iter: Vec<Area>,
    /// again, dz = sin theta dx
    /// so dz/dx = sin theta
    theta_vec: Vec<Angle>,
    /// this is length of each component
    dx: Vec<Length>,
    /// pressure vector 
    pressure_vector_last_iter: Vec<Pressure>,
    /// this is the hydraulic diameter of each component 
    /// used to compute darcy friction factor
    d_h: Vec<Length>,
}
