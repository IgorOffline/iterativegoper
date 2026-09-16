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
use constants::{BACKGROUND, MATERIAL_PALETTE_V3_SHADE_500, SEEDS};
use geometry::Point;
use simulation::Simulation;

fn window_conf() -> Conf {
    Conf {
        window_title: "bravo_leaf".to_owned(),
        window_width: 500,
        window_height: 500,
        ..Default::default()
    }
}

fn is_demo() -> bool {
    std::env::args().skip(1).any(|arg| arg == "demo")
}

async fn run_demo() {
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

async fn run_static() {
    let mut subject = Simulation::new(
        SEEDS[0],
        MATERIAL_PALETTE_V3_SHADE_500[1],
        Point::new(0.0, 0.0),
        (0.0, 0.0),
    );
    subject.grow_to_maturity();
    let mut screenshot_index = 0usize;
    loop {
        clear_background(BACKGROUND);
        subject.draw();
        if is_key_pressed(KeyCode::P) {
            let path = format!("screenshot_{screenshot_index:04}.png");
            get_screen_data().export_png(&path);
            screenshot_index += 1;
        }
        next_frame().await;
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    if is_demo() {
        run_demo().await;
    } else {
        run_static().await;
    }
}
