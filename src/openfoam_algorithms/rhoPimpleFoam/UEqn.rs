use tuas_boussinesq_solver::fluid_mechanics_correlations::pipe_calculations::pipe_calc_pressure_loss;
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
    /// this is wetted perimeter of each component
    p_w: Vec<Length>,
    /// pressure vector 
    pressure_vector_last_iter: Vec<Pressure>,
    /// this is the hydraulic diameter of each component 
    /// used to compute darcy friction factor
    d_h: Vec<Length>,
    /// this is the fluid viscosity
    mu: Vec<DynamicViscosity>,
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
    /// In this case, A is in units of frequency
    ///
    /// upwinding is used to interpolate flux at the faces
    /// note: this assumes the control volumes are connected in a loop
    pub fn A(&self,
        include_transient_term: bool,
        include_advection_term: bool,
    ) -> Vec<Frequency> {
        let vector_length: usize = self.mass_flowrate_vector_last_iter.len();
        let mut a_vec: Vec<Frequency> 
            = vec![Frequency::ZERO; vector_length];
        // obtain timestep
        let delta_t = self.delta_t;

        for (i, mass_rate_ptr) in self.mass_flowrate_vector_last_iter.iter().enumerate(){


            let rho: MassDensity = self.density_vector_last_iter[i];
            // these are not needed, just commenting out
            //let theta: Angle = self.theta_vec[i];
            //let g: Acceleration = self.g;
            let xs_area: Area = self.cross_sectional_area_vec_last_iter[i];

            let hydraulic_diameter = self.d_h[i];
            let pipe_length = self.dx[i];

            // let's consider the fanning term: 
            //
            // delta P = 1/2 rho u^2 (f L/D + K) 
            //
            // delta P = 1/2 rho (m/(rho A_XS))^2 (f L/D + K)
            // delta P = 1/2 rho (m^2/(rho^2 A_XS^2)) (f L/D + K)
            // delta P = 1/2 (m^2/(rho A_XS^2)) (f L/D + K)
            //
            // we multiply by wetted perimeter
            // delta P * p_w = 1/2 (m^2/(rho A_XS^2)) * p_w (f L/D + K)
            //
            // Take it on a per unit mass basis
            //
            // delta P * p_w/m = 1/2 (m/(rho A_XS^2)) * p_w (f L/D + K)
            //
            // we then get our answer for getting the term: 
            //
            // 1/2 (m/(rho A_XS^2)) * p_w (f L/D + K)

            let absolute_mass_flowrate: MassRate = mass_rate_ptr.abs();
            let fluid_viscosity = self.mu[i];
            let wetted_perimeter: Length = self.p_w[i];

            // note that for steam turbine, we just assume zero K term 
            // and smooth pipe, for ease of calculation

            let form_loss_k: Ratio = Ratio::ZERO;
            let absolute_roughness: Length = Length::ZERO;

            let pressure_drop: Pressure = 
                pipe_calc_pressure_loss(
                    absolute_mass_flowrate, 
                    xs_area, 
                    hydraulic_diameter, 
                    fluid_viscosity, 
                    rho, 
                    pipe_length, 
                    absolute_roughness, 
                    form_loss_k.into()
                ).unwrap();

            // this is the friction loss term
            let mut A: Frequency = 
                pressure_drop/
                absolute_mass_flowrate*
                wetted_perimeter;


            // now, if we don't include transient or advection term, 
            // just push A to a_vector 

            if !include_transient_term && !include_advection_term {
                a_vec[i] = A;
                continue;
            };

            // now, let's see if we include the transient term, 
            // this will be quite simple 

            if include_transient_term {
                A += delta_t.recip();
            };

            // if we don't include the advection term, then exit 

            if !include_advection_term {
                a_vec[i] = A;
                continue;
            };

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

            let advection_term: Frequency = 
                total_outflow_from_cell/
                pipe_length/ 
                xs_area/
                rho;


            // we add the advection term in
            A += advection_term;

            
            // then we return A
            a_vec[i] = A;

        }

        // after all this, we return the a_vec 

        return a_vec;


    }

    /// this is pressure gradient term 
    /// A_XS * (0.5 P_CV+1 - 0.5 P_CV-1)/delta_x
    pub fn get_dpdx_term_times_area_times_delta_t(&self) -> Vec<MassRate>{

        // obtain timestep
        let delta_t = self.delta_t;

        let vector_length: usize = self.mass_flowrate_vector_last_iter.len();
        let mut dpdx_term_times_area_times_delta_t: Vec<MassRate> = 
            vec![MassRate::ZERO; vector_length];
        let final_index: usize = self.mass_flowrate_vector_last_iter.len() - 1;
        

        for (i, dx_ptr) in self.dx.iter().enumerate(){

            let dx: Length = *dx_ptr;
            let xs_area: Area = self.cross_sectional_area_vec_last_iter[i];

            let p_cv_plus_one: Pressure;
            let p_cv_minus_one: Pressure;

            if i == 0 {
                // if this is the first control vol in the vector
                // we need to handle the edge case
                // and go like in a cyclic BC
                // This is most definitely true for Rankine Cycle
                p_cv_minus_one = self.pressure_vector_last_iter[final_index];
                p_cv_plus_one = self.pressure_vector_last_iter[i+1];
            } else if i == final_index {
                // if this is the last control volume in the vector, we also need 
                // to handle that edge case, 
                // this assumes that control volumes are in one straight line in series
                //
                // and go like in a cyclic BC
                // This is most definitely true for Rankine Cycle
                p_cv_minus_one = self.pressure_vector_last_iter[i-1];
                p_cv_plus_one = self.pressure_vector_last_iter[0];
            } else {
                p_cv_minus_one = self.pressure_vector_last_iter[i-1];
                p_cv_plus_one = self.pressure_vector_last_iter[i+1];
            }

            let dp: Pressure = 0.5 * p_cv_plus_one - 0.5 * p_cv_minus_one;
            let dpdx = dp/dx;


            dpdx_term_times_area_times_delta_t[i] = 
                dpdx * xs_area * delta_t;

        }

        return dpdx_term_times_area_times_delta_t;
    }
    


}
