use macroquad::prelude::*;

mod blade;
mod constants;
mod geometry;
mod growth_model;
mod rendering;
mod simulation;
mod venation;

use blade::random_blade_point;
use constants::{MATERIAL_PALETTE_V3_SHADE_500, SEEDS};
use geometry::Point;
use simulation::Simulation;

fn roll_seed() -> u64 {
    let roll = ::rand::random_range(0..SEEDS.len());
    SEEDS[roll]
}

fn roll_color() -> Color {
    let roll = ::rand::random_range(0..MATERIAL_PALETTE_V3_SHADE_500.len());
    MATERIAL_PALETTE_V3_SHADE_500[roll]
}

fn roll_root() -> Point {
    random_blade_point()
}

fn roll_screen_offset() -> (f32, f32) {
    (
        ::rand::random_range(-175.0..=175.0),
        ::rand::random_range(-175.0..=175.0),
    )
}

fn window_conf() -> Conf {
    Conf {
        window_title: "bravo_leaf".to_owned(),
        window_width: 500,
        window_height: 500,
        ..Default::default()
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let mut simulation = Simulation::new(
        SEEDS[0],
        MATERIAL_PALETTE_V3_SHADE_500[2],
        roll_root(),
        roll_screen_offset(),
    );
    loop {
        if is_key_pressed(KeyCode::R) {
            simulation =
                Simulation::new(roll_seed(), roll_color(), roll_root(), roll_screen_offset());
        }
        simulation.step(get_frame_time() as f64);
        simulation.draw();
        next_frame().await;
    }
}
