use uom::ConstZero;
use uom::si::angle::degree;
use uom::si::area::square_meter;
use uom::si::energy::joule;
use uom::si::magnetic_flux_density::tesla;
use uom::si::moment_of_inertia::kilogram_square_meter;
use uom::si::ratio::ratio;
use uom::si::torque::newton_meter;
/// This is a 3 phase generator where 
/// phase shift is 0 degrees, 60 degrees and 120 degrees 
/// respectively
#[allow(non_snake_case)]
use uom::si::f64::*;
#[derive(Debug,Clone, PartialEq)]
pub struct ThreePhaseElectricGeneratorTurbine {
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

/// these are defaults for a three phase generator
impl ThreePhaseElectricGeneratorTurbine {

    pub fn new_250_megawatt_generator() -> Self {

        let B = MagneticFluxDensity::new::<tesla>(1.0);
        let A = Area::new::<square_meter>(1.0);
        let N: usize = 25;
        let I = MomentOfInertia::new::<kilogram_square_meter>(100_000.0);
        let eta = Ratio::new::<ratio>(0.98);
        let omega = AngularVelocity::ZERO;

        return Self {
            I,
            omega,
            N,
            B,
            A,
            eta,
        };
    }
    pub fn new(
        B: MagneticFluxDensity,
        A: Area,
        N: usize,
        I: MomentOfInertia,
        eta: Ratio,
        omega: AngularVelocity,
    ) -> Self {


        return Self {
            I,
            omega,
            N,
            B,
            A,
            eta,
        };
    }
}

impl ThreePhaseElectricGeneratorTurbine {

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
    ///
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
    /// this mutably calculates new angular velocity 
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
    pub fn advance_timestep(
        &mut self,
        torque_source: Torque,
        load_resistance: ElectricalResistance,
        current_time: Time,
        delta_t: Time,
    ){
        let new_angular_velocity = self.calculate_new_angular_velocity(
            torque_source, load_resistance, current_time, delta_t);

        self.omega = new_angular_velocity;
    }

    pub fn set_magnetic_field(&mut self,B: MagneticFluxDensity){
        self.B = B
    }

    pub fn get_emf_1(&self, t: Time) -> ElectricPotential {

        let nba: MagneticFlux = self.N as f64 * self.B * self.A;
        let omega = self.omega;

        let phase_shift_1 = Angle::ZERO;

        let theta: Angle = (self.omega * t).into();
        let cos_angle_1: Ratio = (theta + phase_shift_1).cos();

        let emf = -nba * omega * cos_angle_1;

        return emf;

    }

    pub fn get_emf_2(&self, t: Time) -> ElectricPotential {

        let nba: MagneticFlux = self.N as f64 * self.B * self.A;
        let omega = self.omega;

        let phase_shift_2 = Angle::new::<degree>(120.0);

        let theta: Angle = (self.omega * t).into();
        let cos_angle_2: Ratio = (theta + phase_shift_2).cos();

        let emf = -nba * omega * cos_angle_2;

        return emf;

    }
    pub fn get_emf_3(&self, t: Time) -> ElectricPotential {

        let nba: MagneticFlux = self.N as f64 * self.B * self.A;
        let omega = self.omega;

        let phase_shift_3 = Angle::new::<degree>(240.0);

        let theta: Angle = (self.omega * t).into();
        let cos_angle_3: Ratio = (theta + phase_shift_3).cos();

        let emf = -nba * omega * cos_angle_3;

        return emf;

    }

    pub fn get_power(&self,
        load_resistance: ElectricalResistance,
        t: Time
    ) -> Power {

        let emf_1 = self.get_emf_1(t);
        let emf_2 = self.get_emf_2(t);
        let emf_3 = self.get_emf_3(t);

        let p: Power = load_resistance.recip() *
            (emf_1 * emf_1
            +emf_2 * emf_2
            +emf_3 * emf_3);

        return p;
            
    }

    pub fn get_current_1(&self, 
        load_resistance: ElectricalResistance,
        t: Time) -> ElectricCurrent {
        
        self.get_emf_1(t)/load_resistance
    }

    pub fn get_current_2(&self, 
        load_resistance: ElectricalResistance,
        t: Time) -> ElectricCurrent {
        
        self.get_emf_2(t)/load_resistance
    }
    pub fn get_current_3(&self, 
        load_resistance: ElectricalResistance,
        t: Time) -> ElectricCurrent {
        
        self.get_emf_3(t)/load_resistance
    }

    /// sets the rpm
    pub fn set_omega(&mut self, omega: AngularVelocity){
        self.omega = omega
    }

    /// gets the rpm
    pub fn get_omega(&self) -> AngularVelocity{
        self.omega
    }

}
