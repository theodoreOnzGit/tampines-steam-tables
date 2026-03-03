use egui::Color32;
/// From ChatGPT
/// Steps:
/// Cold colors (blue) start with high values in the blue channel (B = 1, G = 0).
/// Hot colors (red) end with high values in the red channel (R = 1, G = 0).
pub fn hot_to_cold_colour_mark_1(hotness: f32) -> Color32 {
    let mut hotness_clone = hotness.clone();

    // ensures hotness is between 0 and 1
    if hotness_clone < 0.0 {
        hotness_clone = 0.0;
    } else if hotness_clone > 1.0 {
        hotness_clone = 1.0
    }

    let red: f32 = 255.0 * hotness_clone;
    let green: f32 = 135.0 * (1.0 - hotness_clone);
    let blue: f32 = 255.0 * (1.0 - hotness_clone);

    return Color32::from_rgb(
        red as u8, 
        green as u8, 
        blue as u8);
}

/// from colour picker 
/// from black(cold) to red (hot)
pub fn hot_to_cold_colour_mark_2(hotness: f32) -> Color32 {
    let mut hotness_clone = hotness.clone();

    // ensures hotness is between 0 and 1
    if hotness_clone < 0.0 {
        hotness_clone = 0.0;
    } else if hotness_clone > 1.0 {
        hotness_clone = 1.0
    }

    let red: f32 = 255.0 * hotness_clone;
    let green: f32 = 0.0;
    let blue: f32 = 0.0;

    return Color32::from_rgb(
        red as u8, 
        green as u8, 
        blue as u8);
}

