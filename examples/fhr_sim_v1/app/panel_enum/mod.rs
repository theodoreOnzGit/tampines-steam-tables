
#[derive(serde::Deserialize, serde::Serialize,PartialEq,Clone, Debug)]
pub enum Panel {
    MainPage,
    ReactorPowerGraphs,
    PoisonGraphs,

}
