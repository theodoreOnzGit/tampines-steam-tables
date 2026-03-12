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
    /// upwinding is used to interpolate flux at the faces
    pub fn A(&self,
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



            // now, if we consider the advection and transient term,
            // we need to get the face flux 

            // case 1: positive flow 1 -> 2 -> 3
            // or 
            // (n-1) -> (n) -> (n+1)

            let m_n_flow: MassRate = *mass_rate_ptr;

            // now, for m_n_minus_one and m_n_plus_one flow, we need 
            // to handle edge cases

            let m_n_minus_one_flow: MassRate;
            let m_n_plus_one_flow: MassRate;
            let dx_n: Length = self.dx[i];
            let dx_n_minus_one: Length;
            let dx_n_plus_one: Length;
            let final_index: usize = self.mass_flowrate_vector_last_iter.len() - 1;

            // if this is the first control vol in the vector
            // we need to handle the edge case
            // and go like in a cyclic BC
            // This is most definitely true for Rankine Cycle
            if i == 0 {
                m_n_minus_one_flow = self.mass_flowrate_vector_last_iter[final_index];
                dx_n_minus_one = self.dx[final_index];
            } else {
                m_n_minus_one_flow = self.mass_flowrate_vector_last_iter[i-1];
                dx_n_minus_one = self.dx[i-1];
            }

            // if this is the last control volume in the vector, we also need 
            // to handle that edge case, 
            // this assumes that control volumes are in one straight line in series
            //
            // and go like in a cyclic BC
            // This is most definitely true for Rankine Cycle
            if i == final_index {
                m_n_plus_one_flow = self.mass_flowrate_vector_last_iter[0];
                dx_n_plus_one = self.dx[0];
            } else {
                m_n_plus_one_flow = self.mass_flowrate_vector_last_iter[i+1];
                dx_n_plus_one = self.dx[i+1];
            }


            // now let's handle case 1, 2 and 3
            //
            let m_n_minus_one_flow_positive: bool = 
                m_n_minus_one_flow > MassRate::ZERO;
            let m_n_flow_positive: bool = 
                m_n_flow > MassRate::ZERO;
            let m_n_plus_one_flow_positive: bool = 
                m_n_plus_one_flow > MassRate::ZERO;

            let m_n_minus_1_n_face_flow: MassRate;
            let m_n_n_plus_1_face_flow: MassRate;
            let total_length_n_minus_1_n: Length = dx_n + dx_n_minus_one;
            let total_length_n_n_plus_1: Length = dx_n + dx_n_plus_one;
            

            // let's do the face at (n-1) | (n)

            if m_n_minus_one_flow_positive && m_n_flow_positive {
                m_n_minus_1_n_face_flow = m_n_minus_one_flow;
            } else if !m_n_minus_one_flow_positive && !m_n_flow_positive {
                m_n_minus_1_n_face_flow = m_n_flow;
            } else {
                m_n_minus_1_n_face_flow 
                    = m_n_minus_one_flow * dx_n/total_length_n_minus_1_n 
                    + m_n_flow * dx_n_minus_one/total_length_n_minus_1_n;
            }

            // let's do the face at (n) | (n+1)
            if m_n_plus_one_flow_positive && m_n_flow_positive {
                m_n_n_plus_1_face_flow = m_n_flow;
            } else if !m_n_plus_one_flow_positive && !m_n_flow_positive {
                m_n_n_plus_1_face_flow = m_n_flow;
            } else {
                m_n_n_plus_1_face_flow 
                    = m_n_plus_one_flow * dx_n/total_length_n_n_plus_1 
                    + m_n_flow * dx_n_plus_one/total_length_n_n_plus_1;
            }

            // now we can compute total flow out from cell,
            // which is the meaning of the divergence operator

            let total_outflow_from_cell: MassRate = 
                -m_n_minus_1_n_face_flow
                +m_n_n_plus_1_face_flow;

        }


    }



}

