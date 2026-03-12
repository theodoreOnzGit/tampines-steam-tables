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
use uom::{ConstZero, si::{angle::radian, f64::*}};

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
    ///
    /// theta is the angle of incline
    theta_vec: Vec<Angle>,
    /// this is length of each component
    dx: Vec<Length>,
    /// pressure vector 
    pressure_vector_last_iter: Vec<Pressure>,
    /// this is the hydraulic diameter of each component 
    /// used to compute darcy friction factor
    d_h: Vec<Length>,
}


impl UEqn {
    // these will be constructors

}


impl UEqn {

    /// this gives a vector of the velocities at the next iteration 
    /// based on the supplied parameters
    pub fn U(&self) -> Vec<Velocity>{

        todo!()
    }


    /// this gives a vector of H
    ///
    /// H = mass_flowrate/delta_t 
    /// - rho g A_{XS} dz/dx
    ///
    /// there is no common unit for mass flowrate/delta_t 
    /// so i rather return 
    ///
    /// H * delta_t
    pub fn H_times_delta_t(&self, include_transient_term: bool) -> Vec<MassRate> {

        let vector_length: usize = self.mass_flowrate_vector_last_iter.len();
        let mut h_times_delta_t_vec: Vec<MassRate> 
            = vec![MassRate::ZERO; vector_length];
        // obtain timestep
        let delta_t = self.delta_t;
        for (i, mass_rate_ptr) in self.mass_flowrate_vector_last_iter.iter().enumerate(){

            let theta: Angle = self.theta_vec[i];

            let rho: MassDensity = self.density_vector_last_iter[i];
            let g: Acceleration = self.g;
            let xs_area: Area = self.cross_sectional_area_vec_last_iter[i];

            // dz/dx = sin theta, where theta is the angle of incline
            let dz_dx: Ratio = theta.get::<radian>().sin().into();

            let mut h = - rho * g * xs_area * dz_dx;

            
            if include_transient_term {

                h += *mass_rate_ptr/delta_t;

            }

            h_times_delta_t_vec[i] = h * delta_t;



        }
        return h_times_delta_t_vec;
    }



    /// returns the A term, using an upwind scheme to interpolate mass 
    /// flowrates at the faces
    ///
    ///
    /// upwinding scenarios are described in documents
    pub fn A_upwind(&self,
        include_transient_term: bool,
        include_advection_term: bool,
    ){
        let vector_length: usize = self.mass_flowrate_vector_last_iter.len();
        let mut a_vec: Vec<Frequency> 
            = vec![Frequency::ZERO; vector_length];
        // obtain timestep
        let delta_t = self.delta_t;

        for (i, mass_rate_ptr) in self.mass_flowrate_vector_last_iter.iter().enumerate(){

            let theta: Angle = self.theta_vec[i];

            let rho: MassDensity = self.density_vector_last_iter[i];
            let g: Acceleration = self.g;
            let xs_area: Area = self.cross_sectional_area_vec_last_iter[i];


        }


    }

}

