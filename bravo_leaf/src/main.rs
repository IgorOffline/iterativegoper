use macroquad::prelude::*;

mod blade;
mod constants;
mod geometry;
mod growth_model;
mod rendering;
mod simulation;
mod venation;

use constants::{MATERIAL_PALETTE_V3_SHADE_500, SEEDS};
use simulation::Simulation;

fn roll_seed() -> u64 {
    let roll = ::rand::random_range(0..SEEDS.len());
    SEEDS[roll]
}

fn roll_color() -> Color {
    let roll = ::rand::random_range(0..MATERIAL_PALETTE_V3_SHADE_500.len());
    MATERIAL_PALETTE_V3_SHADE_500[roll]
}

#[macroquad::main("bravo_leaf")]
async fn main() {
    let mut simulation = Simulation::new(SEEDS[0], MATERIAL_PALETTE_V3_SHADE_500[2]);
    loop {
        if is_key_pressed(KeyCode::R) {
            simulation = Simulation::new(roll_seed(), roll_color());
        }
        simulation.step(get_frame_time() as f64);
        simulation.draw();
        next_frame().await;
    }
}
