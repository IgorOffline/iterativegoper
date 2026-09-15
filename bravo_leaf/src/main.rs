use macroquad::prelude::*;

mod blade;
mod constants;
mod geometry;
mod growth_model;
mod rendering;
mod simulation;
mod venation;

use simulation::Simulation;

#[macroquad::main("bravo_leaf")]
async fn main() {
    let mut simulation = Simulation::new();
    loop {
        if is_key_pressed(KeyCode::R) {
            simulation = Simulation::new();
        }
        simulation.step(get_frame_time() as f64);
        simulation.draw();
        next_frame().await;
    }
}
