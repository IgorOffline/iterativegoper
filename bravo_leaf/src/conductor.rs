use macroquad::prelude::clear_background;

use crate::blade::random_blade_point;
use crate::constants::*;
use crate::simulation::Simulation;

pub struct Conductor {
    subjects: Vec<Simulation>,
    time_since_spawn: f64,
}

impl Conductor {
    pub fn new() -> Self {
        Conductor {
            subjects: vec![spawn_subject()],
            time_since_spawn: 0.0,
        }
    }

    pub fn step(&mut self, delta_time: f64) {
        self.time_since_spawn += delta_time;
        while self.time_since_spawn >= SUBJECT_SPAWN_INTERVAL {
            self.time_since_spawn -= SUBJECT_SPAWN_INTERVAL;
            self.subjects.push(spawn_subject());
        }

        for subject in self.subjects.iter_mut() {
            subject.step(delta_time);
        }
        self.subjects.retain(|subject| !subject.is_expired());
    }

    pub fn draw(&self) {
        clear_background(BACKGROUND);
        for subject in &self.subjects {
            subject.draw();
        }
    }
}

fn spawn_subject() -> Simulation {
    let seed = SEEDS[::rand::random_range(0..SEEDS.len())];
    let color =
        MATERIAL_PALETTE_V3_SHADE_500[::rand::random_range(0..MATERIAL_PALETTE_V3_SHADE_500.len())];
    let screen_offset = (
        ::rand::random_range(-175.0..=175.0),
        ::rand::random_range(-175.0..=175.0),
    );
    Simulation::new(seed, color, random_blade_point(), screen_offset)
}
