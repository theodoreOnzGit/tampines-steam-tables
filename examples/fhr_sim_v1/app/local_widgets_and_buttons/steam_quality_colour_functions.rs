use egui::Color32;

/// steam quality colour from 0 to 1 
/// blue is zero 
/// white is 1
pub fn steam_quality_colour_mark_1(steam_quality: f32) -> Color32 {
    let mut steam_quality_clone = steam_quality.clone();

    // ensures hotness is between 0 and 1
    if steam_quality_clone < 0.0 {
        steam_quality_clone = 0.0;
    } else if steam_quality_clone > 1.0 {
        steam_quality_clone = 1.0
    }

    let blue: f32 = 255.0 * steam_quality_clone;
    let green: f32 = 255.0 * (1.0 - steam_quality_clone);
    let red: f32 = 255.0 * (1.0 - steam_quality_clone);

    return Color32::from_rgb(
        red as u8, 
        green as u8, 
        blue as u8);
}

