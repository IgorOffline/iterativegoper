use macroquad::prelude::{Color, draw_line, rand};

use crate::constants::*;
use crate::geometry::Point;
use crate::growth_model::{GrowthTransition, VerticalMap, propagate, reference_bounding_box};
use crate::rendering::{fit_world, to_screen};
use crate::venation::{
    RadiusScale, VeinGraph, murray_radii, throw_darts, throw_darts_current, venation_step,
};

pub struct Simulation {
    graph: VeinGraph,
    sources: Vec<Point>,
    time: f64,
    accumulated_time: f64,
    color: Color,
    screen_offset: (f32, f32),
    lifetime: f64,
}

impl Simulation {
    pub fn new(seed: u64, color: Color, root: Point, screen_offset: (f32, f32)) -> Self {
        rand::srand(seed);
        let graph = VeinGraph::seed(root);
        let sources = throw_darts(&graph.nodes, DART_ATTEMPTS);
        Simulation {
            graph,
            sources,
            time: T0,
            accumulated_time: 0.0,
            color,
            screen_offset,
            lifetime: 0.0,
        }
    }

    fn grow_cycle(&mut self) {
        if self.time >= T_MAX {
            return;
        }
        let time_from = self.time;
        let time_to = (time_from + DEV_DT).min(T_MAX);
        let map_from = VerticalMap::build(time_from);
        let map_to = VerticalMap::build(time_to);
        let transition = GrowthTransition {
            time_from,
            time_to,
            map_from: &map_from,
            map_to: &map_to,
        };

        for node in self.graph.nodes.iter_mut() {
            *node = propagate(*node, &transition);
        }
        for source in self.sources.iter_mut() {
            *source = propagate(*source, &transition);
        }
        self.time = time_to;

        let fresh = throw_darts_current(
            &self.graph,
            &self.sources,
            time_to,
            &map_to,
            reference_bounding_box(time_to),
        );
        self.sources.extend(fresh);

        for _ in 0..VENATION_SUBSTEPS {
            venation_step(&mut self.graph, &mut self.sources);
        }
    }

    pub fn step(&mut self, delta_time: f64) {
        self.lifetime += delta_time;
        if self.time >= T_MAX {
            return;
        }
        self.accumulated_time += delta_time;
        while self.accumulated_time >= AUTO_STEP_INTERVAL {
            self.accumulated_time -= AUTO_STEP_INTERVAL;
            self.grow_cycle();
        }
    }

    pub fn is_expired(&self) -> bool {
        self.lifetime >= SUBJECT_LIFESPAN
    }

    fn fade_alpha(&self) -> f32 {
        let remaining = SUBJECT_LIFESPAN - self.lifetime;
        (remaining / SUBJECT_FADE).clamp(0.0, 1.0) as f32
    }

    pub fn draw(&self) {
        let mut color = self.color;
        color.a = self.fade_alpha();
        let radii = murray_radii(&self.graph);
        let scale = RadiusScale::new(radii.iter().cloned().fold(0.0f64, f64::max));
        let fit = fit_world(reference_bounding_box(T_MAX))
            .translated(self.screen_offset.0, self.screen_offset.1 - 105.0);

        for &(from, to) in &self.graph.edges {
            let radius = radii[to];
            if !scale.is_secondary(radius) {
                continue;
            }
            let (from_x, from_y) = to_screen(self.graph.nodes[from], &fit);
            let (to_x, to_y) = to_screen(self.graph.nodes[to], &fit);
            draw_line(
                from_x,
                from_y,
                to_x,
                to_y,
                (radius as f32 * VEIN_WIDTH_PX).max(1.0),
                color,
            );
        }
    }
}
