use macroquad::prelude::{screen_height, screen_width};

use crate::geometry::{BoundingBox, Point};

pub struct ViewportFit {
    scale: f32,
    origin_x: f32,
    origin_y: f32,
}

pub fn fit_world(bounds: BoundingBox) -> ViewportFit {
    let margin = 60.0f32;
    let available_width = screen_width() - 2.0 * margin;
    let available_height = screen_height() - 2.0 * margin;
    let span_x = (bounds.max_x - bounds.min_x).max(1e-6) as f32;
    let span_y = (bounds.max_y - bounds.min_y).max(1e-6) as f32;
    let scale = (available_width / span_x).min(available_height / span_y);
    let center_x = 0.5 * (bounds.min_x + bounds.max_x) as f32;
    ViewportFit {
        scale,
        origin_x: screen_width() * 0.5 - center_x * scale,
        origin_y: screen_height() - margin,
    }
}

pub fn to_screen(point: Point, fit: &ViewportFit) -> (f32, f32) {
    (
        fit.origin_x + point.x as f32 * fit.scale,
        fit.origin_y - point.y as f32 * fit.scale,
    )
}
