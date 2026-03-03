use uom::si::f64::*;
impl super::TampinesSteamTableCV {

    pub fn set_tpx(&mut self,
        t: ThermodynamicTemperature,
        p: Pressure,
        x: f64){

        let volume = self.get_volume();

        *self = Self::new_from_tp_quality(t, p, volume, x);
        

    }

    pub fn set_ph(&mut self,
        p: Pressure,
        h: AvailableEnergy,){

        let volume = self.get_volume();

        *self = Self::new_from_ph(p, h, volume);
        

    }
    pub fn set_ps(&mut self,
        p: Pressure,
        s: SpecificHeatCapacity,){

        let volume = self.get_volume();

        *self = Self::new_from_ps(p, s, volume);
        

    }
}
