use macroquad::prelude::*;

mod blade;
mod conductor;
mod constants;
mod geometry;
mod growth_model;
mod rendering;
mod simulation;
mod venation;

use conductor::Conductor;

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
    let mut conductor = Conductor::new();
    loop {
        if is_key_pressed(KeyCode::R) {
            conductor = Conductor::new();
        }
        conductor.step(get_frame_time() as f64);
        conductor.draw();
        next_frame().await;
    }
}
