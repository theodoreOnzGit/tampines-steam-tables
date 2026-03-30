use uom::si::f64::*;
use uom::si::force::newton;
use uom::si::pressure::{bar, pascal};
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::ps_flash_eqm::v_ps_eqm;
use crate::prelude::{TampinesSteamTableCV};

#[inline]
pub fn force_balance_isentropic_nozzle(
    p1: Pressure,
    p2: Pressure,
    h1: AvailableEnergy,
    mass_flowrate: MassRate,
    a1: Area,
    a2: Area,
    ) -> Force {
    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1: TampinesSteamTableCV = 
        TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);

    let rho1: MassDensity = state_1.get_specific_volume().recip();
    let s1: SpecificHeatCapacity = state_1.get_specific_entropy();
    // note: this is isentropic, hence s2 = s1;
    let s2 = s1;
    let v1: Velocity = mass_flowrate/rho1/a1;

    let p1a1: Force = p1*a1;
    // left hand side of momentum balance
    //
    // P1 A1 + dot{m}^2/{rho1 a1}
    let lhs: Force = p1a1 + mass_flowrate * v1;

    // if we have p2, we can get the second thermodynamic state
    let p2a2 = p2 * a2;
    let rho2 = v_ps_eqm(p2, s2).recip();
    let rhs: Force = p2a2 + mass_flowrate * mass_flowrate/rho2/a2;

    return lhs-rhs;

}

#[inline]
pub fn print_graph_pts_for_outlet_pressure_and_force_balance(
    p1: Pressure,
    h1: AvailableEnergy,
    a1: Area,
    a2: Area,
    mass_flowrate: MassRate,
){
    // now let's plot a trend of isentropic nozzles 
    println!("{:?}",&("oulet pressure (Bar)","force bal(Newton)"));
    let mut p2 = p1;
    let n = 100;
    for i in 0..n {

        p2 = (n as f64 - i as f64)/(n as f64) * p1;

        let force_bal: Force = 
            force_balance_isentropic_nozzle(p1, p2, h1, mass_flowrate, a1, a2);

        println!("{:?}",&(p2.get::<bar>(),force_bal.get::<newton>()));

    }

    // let's also check the lower bound pressures
    // below the first 100 points
    for i in 0..n {

        if p2 > Pressure::new::<pascal>(1000_f64) {

            p2 = (1_f64)/(n as f64) * p1;
            p2 *= (n-i) as f64 / (n as f64);


            let force_bal: Force = 
                force_balance_isentropic_nozzle(p1, p2, h1, mass_flowrate, a1, a2);

            println!("{:?}",&(p2.get::<bar>(),force_bal.get::<newton>()));

        }

    }
}

