use macroquad::rand;

use crate::blade::{leaf_half_width, point_in_blade};
use crate::constants::*;
use crate::geometry::{BoundingBox, Point, distance_squared};
use crate::growth_model::{VerticalMap, unpropagate};

pub struct VeinGraph {
    pub nodes: Vec<Point>,
    pub edges: Vec<(usize, usize)>,
}

impl VeinGraph {
    pub fn seed(root: Point) -> Self {
        VeinGraph {
            nodes: vec![root],
            edges: Vec::new(),
        }
    }
}

fn nearest_node(nodes: &[Point], target: Point) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (index, &node) in nodes.iter().enumerate() {
        let distance = distance_squared(node, target);
        if best.map_or(true, |(_, best_distance)| distance < best_distance) {
            best = Some((index, distance));
        }
    }
    best.map(|(index, _)| index)
}

pub fn throw_darts(nodes: &[Point], attempts: usize) -> Vec<Point> {
    let source_spacing_squared = BIRTH_DIST_SOURCE * BIRTH_DIST_SOURCE;
    let vein_spacing_squared = BIRTH_DIST_VEIN * BIRTH_DIST_VEIN;
    let mut sources: Vec<Point> = Vec::new();
    for _ in 0..attempts {
        let height = rand::gen_range(0.0, Y0_MAX);
        let half_width = leaf_half_width(height);
        let candidate = Point::new(rand::gen_range(-half_width, half_width), height);
        if !point_in_blade(candidate) {
            continue;
        }
        if sources
            .iter()
            .any(|&source| distance_squared(source, candidate) < source_spacing_squared)
        {
            continue;
        }
        if nodes
            .iter()
            .any(|&node| distance_squared(node, candidate) < vein_spacing_squared)
        {
            continue;
        }
        sources.push(candidate);
    }
    sources
}

pub fn throw_darts_current(
    graph: &VeinGraph,
    sources: &[Point],
    time: f64,
    map: &VerticalMap,
    bounds: BoundingBox,
) -> Vec<Point> {
    let source_spacing_squared = BIRTH_DIST_SOURCE * BIRTH_DIST_SOURCE;
    let vein_spacing_squared = BIRTH_DIST_VEIN * BIRTH_DIST_VEIN;
    let mut fresh: Vec<Point> = Vec::new();
    for _ in 0..DART_ATTEMPTS_CYCLE {
        let candidate = Point::new(
            rand::gen_range(bounds.min_x, bounds.max_x),
            rand::gen_range(bounds.min_y, bounds.max_y),
        );
        if !point_in_blade(unpropagate(candidate, time, map)) {
            continue;
        }
        if sources
            .iter()
            .any(|&source| distance_squared(source, candidate) < source_spacing_squared)
        {
            continue;
        }
        if fresh
            .iter()
            .any(|&source| distance_squared(source, candidate) < source_spacing_squared)
        {
            continue;
        }
        if graph
            .nodes
            .iter()
            .any(|&node| distance_squared(node, candidate) < vein_spacing_squared)
        {
            continue;
        }
        fresh.push(candidate);
    }
    fresh
}

pub fn venation_step(graph: &mut VeinGraph, sources: &mut Vec<Point>) {
    if sources.is_empty() || graph.nodes.is_empty() {
        return;
    }
    let kill_distance_squared = KILL_DISTANCE * KILL_DISTANCE;

    let mut influence = vec![(0.0f64, 0.0f64); graph.nodes.len()];
    let mut has_influence = vec![false; graph.nodes.len()];
    for &source in sources.iter() {
        if let Some(node_index) = nearest_node(&graph.nodes, source) {
            let node_position = graph.nodes[node_index];
            let (delta_x, delta_y) = (source.x - node_position.x, source.y - node_position.y);
            let length = (delta_x * delta_x + delta_y * delta_y).sqrt();
            if length > 1e-12 {
                influence[node_index].0 += delta_x / length;
                influence[node_index].1 += delta_y / length;
                has_influence[node_index] = true;
            }
        }
    }

    let original_node_count = graph.nodes.len();
    for node_index in 0..original_node_count {
        if !has_influence[node_index] {
            continue;
        }
        let (accumulated_x, accumulated_y) = influence[node_index];
        let magnitude = (accumulated_x * accumulated_x + accumulated_y * accumulated_y).sqrt();
        if magnitude <= 1e-12 {
            continue;
        }
        let node_position = graph.nodes[node_index];
        let new_position = Point::new(
            node_position.x + VEIN_GROWTH_STEP * accumulated_x / magnitude,
            node_position.y + VEIN_GROWTH_STEP * accumulated_y / magnitude,
        );
        let new_index = graph.nodes.len();
        graph.nodes.push(new_position);
        graph.edges.push((node_index, new_index));
    }

    sources.retain(|&source| {
        graph
            .nodes
            .iter()
            .all(|&node| distance_squared(node, source) > kill_distance_squared)
    });
}

pub fn murray_radii(graph: &VeinGraph) -> Vec<f64> {
    let node_count = graph.nodes.len();
    let mut child_power_sum = vec![0.0f64; node_count];
    let mut has_child = vec![false; node_count];
    let mut parent = vec![None::<usize>; node_count];
    for &(from, to) in &graph.edges {
        parent[to] = Some(from);
    }
    let mut radius = vec![TIP_RADIUS; node_count];
    for node_index in (0..node_count).rev() {
        radius[node_index] = if has_child[node_index] {
            child_power_sum[node_index].powf(1.0 / MURRAY_EXPONENT)
        } else {
            TIP_RADIUS
        };
        if let Some(parent_index) = parent[node_index] {
            child_power_sum[parent_index] += radius[node_index].powf(MURRAY_EXPONENT);
            has_child[parent_index] = true;
        }
    }
    radius
}

pub struct RadiusScale {
    maximum: f64,
}

impl RadiusScale {
    pub fn new(maximum_radius: f64) -> Self {
        RadiusScale {
            maximum: maximum_radius,
        }
    }

    pub fn is_secondary(&self, radius: f64) -> bool {
        if self.maximum <= 0.0 {
            return false;
        }
        let fraction = radius / self.maximum;
        (1.0 / 3.0..2.0 / 3.0).contains(&fraction)
    }
}
