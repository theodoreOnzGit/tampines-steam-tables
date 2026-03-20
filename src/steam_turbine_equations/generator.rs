use uom::ConstZero;
use uom::si::angle::degree;
use uom::si::energy::joule;
use uom::si::torque::newton_meter;
/// This is a 3 phase generator where 
/// phase shift is 0 degrees, 60 degrees and 120 degrees 
/// respectively
#[allow(non_snake_case)]
use uom::si::f64::*;
#[derive(Debug,Clone, PartialEq)]
pub struct ThreePhaseElectricGenerator {
    /// I
    I: MomentOfInertia,

    /// omega 
    omega: AngularVelocity,

    /// no of turns for the coil
    N: usize,

    /// magnetic flux density for generator
    B: MagneticFluxDensity,

    /// area of coil 
    A: Area,

    /// turbine efficiency 
    eta: Ratio,


}

impl ThreePhaseElectricGenerator {

    /// this immutably calculates new angular velocity 
    /// in an explicit manner, given a source term
    ///     \begin{equation*}
    /// 	\omega^{t+ \Delta t} \left(
    /// 		\frac{I}{\Delta t} 
    /// 		+ \frac{( N^2 B^2 A^2 )}{\eta R_{load} } \sum_j \cos^2 (\omega^t t + b_j)  
    /// \right)
    /// 	= 
    /// 	I \frac{\omega^{t } }{\Delta t}
    /// 	+ \text{source} 
    /// \end{equation*}
    ///
    /// omega^{t + Delta t} (I/delta_t + (NBA)^2/(eta R_load) 
    /// sum (cos^2 (omega^t t + b_j))
    pub fn calculate_new_angular_velocity(
        &self,
        source: Torque,
        load_resistance: ElectricalResistance,
        current_time: Time,
        delta_t: Time,
    ) -> AngularVelocity {

        let t = current_time;

        let theta: Angle = (self.omega * t).into();

        let phase_shift_1 = Angle::ZERO;
        let phase_shift_2 = Angle::new::<degree>(120.0);
        let phase_shift_3 = Angle::new::<degree>(240.0);

        let cos_angle_1: Ratio = (theta + phase_shift_1).cos();
        let cos_angle_2: Ratio = (theta + phase_shift_2).cos();
        let cos_angle_3: Ratio = (theta + phase_shift_3).cos();

        let cosine_summation: Ratio = 
            cos_angle_1 * cos_angle_1 +
            cos_angle_2 * cos_angle_2 +
            cos_angle_3 * cos_angle_3;

        let nba: MagneticFlux = self.N as f64 * self.B * self.A;

        let coeff = self.I/delta_t 
            + (nba * nba)
            /self.eta 
            /load_resistance
            * cosine_summation;
        
        let mut rhs: Torque = source;

        let momentum_torque: Energy = self.I * self.omega/delta_t;
        // note that both torque and energy have the units of newton 
        // meter
        //
        // Work done (Joule) = F * ds (also newton meter)
        // uom distinguishes both of them though

        rhs += Torque::new::<newton_meter>(
            momentum_torque.get::<joule>()
        );


        return (rhs/coeff).into();
    }

}
