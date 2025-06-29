
use uom::si::f64::*;
use uom::ConstZero;
use uom::si::power::megawatt;
use uom::si::time::second;

#[derive(Debug,Clone)]
pub struct PagePlotData {
    /// the 
    ///
    /// simulation time, 
    /// reactor power (plus decay heat)
    /// reactor power (without decay heat)
    pub reactor_power_plot_data: Vec<(Time,Power,Power)>,




    // recording interval for graphs
    pub graph_data_record_interval_seconds: f64,

    // recording interval for csv 
    pub csv_display_interval_seconds: f64,




}

pub const NUM_DATA_PTS_IN_PLOTS: usize = 4000;

impl PagePlotData {

    /// inserts a data point, most recent being on top 
    pub fn insert_reactor_power_data(&mut self, 
        simulation_time: Time,
        reactor_power_with_decay_heat: Power,
        reactor_power_without_decay_heat: Power){

        // first convert into a tuple,

        let data_tuple = 
            (simulation_time,reactor_power_with_decay_heat,
             reactor_power_without_decay_heat);

        // now insert this into the heater
        // how?
        // map the vectors out first 
        let mut current_heater_data_vec: Vec< (Time,Power,
            Power)>;

        current_heater_data_vec = self.reactor_power_plot_data.iter().map(|&values|{
            values
        }).collect();

        // now, insert the latest data at the top
        current_heater_data_vec.insert(0, data_tuple);

        // take the first NUM_DATA_PTS_IN_PLOTS pieces as a fixed size array 
        // which is basically the array size

        let mut new_array_to_be_put_back: Vec<(Time,Power,
            Power)> =
            vec![ (Time::ZERO, Power::ZERO, 
             Power::ZERO); NUM_DATA_PTS_IN_PLOTS
            ];

        // map the first NUM_DATA_PTS_IN_PLOTS values of the current heater data vec
        
        for n in 0..NUM_DATA_PTS_IN_PLOTS {
            new_array_to_be_put_back[n] = current_heater_data_vec[n];
        }

        self.reactor_power_plot_data = new_array_to_be_put_back;

    }



    /// gets reactor power data over time
    /// time in second, temp in degc
    pub fn get_reactor_power_megawatts_vs_time_secs_vec(&self) -> Vec<[f64;2]> {

        let time_reactor_power_vec: Vec<[f64;2]> = self.reactor_power_plot_data.iter().map(
            |tuple|{
                let (time,reactor_power,_power_less_decay_heat) = *tuple;

                if reactor_power.get::<megawatt>() > 0.0 {
                    [time.get::<second>(), reactor_power.get::<megawatt>()]
                } else {
                    // don't return anything, a default 0.0 will do 
                    // this is the initial condition

                    // doesn't make sense for reactor power to be less than 0
                    [0.0,0.0]
                }

            }
        ).collect();

        return time_reactor_power_vec;
    }

    /// gets reactor power data over time
    /// time in second, temp in degc
    pub fn get_reactor_power_no_decay_heat_megawatts_vs_time_secs_vec(&self) -> Vec<[f64;2]> {

        let time_reactor_power_vec: Vec<[f64;2]> = self.reactor_power_plot_data.iter().map(
            |tuple|{
                let (time,_reactor_power,reactor_power_less_decay_heat) = *tuple;

                if reactor_power_less_decay_heat.get::<megawatt>() > 0.0 {
                    [time.get::<second>(), reactor_power_less_decay_heat.get::<megawatt>()]
                } else {
                    // don't return anything, a default 0.0 will do 
                    // this is the initial condition

                    // doesn't make sense for reactor power to be less than 0
                    [0.0,0.0]
                }

            }
        ).collect();

        return time_reactor_power_vec;
    }




}

impl Default for PagePlotData {
    fn default() -> Self {

        // basically a whole array of dimensioned zeroes
        let reactor_power_plot_data = 
            vec![ (Time::ZERO, Power::ZERO, 
             Power::ZERO,); NUM_DATA_PTS_IN_PLOTS
            ];

        // by default, record every 0.1s
        let graph_data_record_interval_seconds = 0.1;
        let csv_display_interval_seconds = 0.1;


        Self { 
            // first, a blank dataset
            reactor_power_plot_data,
            graph_data_record_interval_seconds,
            csv_display_interval_seconds,

        }
    }
}


pub mod update;
