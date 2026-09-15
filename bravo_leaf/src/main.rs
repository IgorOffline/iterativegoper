use macroquad::prelude::*;

mod blade;
mod constants;
mod geometry;
mod growth_model;
mod rendering;
mod simulation;
mod venation;

use constants::SEEDS;
use simulation::Simulation;

fn roll_seed() -> u64 {
    let roll = ::rand::random_range(0..SEEDS.len());
    SEEDS[roll]
}

#[macroquad::main("bravo_leaf")]
async fn main() {
    let mut simulation = Simulation::new(SEEDS[0]);
    loop {
        if is_key_pressed(KeyCode::R) {
            simulation = Simulation::new(roll_seed());
        }
        simulation.step(get_frame_time() as f64);
        simulation.draw();
        next_frame().await;
    }
}
