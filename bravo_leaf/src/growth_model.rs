use crate::constants::{N_STEPS, T0, Y0_MAX};
use crate::geometry::{BoundingBox, Point};

fn relative_growth_rate_x(material_height: f64) -> f64 {
    1.30 - 0.20 * material_height
}
fn relative_growth_rate_y(height: f64) -> f64 {
    1.25 - 0.15 * height
}

struct InterpolationAbscissae<'a>(&'a [f64]);
struct InterpolationOrdinates<'a>(&'a [f64]);

fn interpolate(
    abscissae: InterpolationAbscissae,
    ordinates: InterpolationOrdinates,
    query: f64,
) -> f64 {
    let (abscissae, ordinates) = (abscissae.0, ordinates.0);
    let count = abscissae.len();
    if query <= abscissae[0] {
        return ordinates[0];
    }
    if query >= abscissae[count - 1] {
        return ordinates[count - 1];
    }
    let (mut low, mut high) = (0usize, count - 1);
    while high - low > 1 {
        let middle = (low + high) / 2;
        if abscissae[middle] <= query {
            low = middle;
        } else {
            high = middle;
        }
    }
    let fraction = (query - abscissae[low]) / (abscissae[high] - abscissae[low]);
    ordinates[low] + fraction * (ordinates[high] - ordinates[low])
}

pub struct VerticalMap {
    material_grid: Vec<f64>,
    heights: Vec<f64>,
}

impl VerticalMap {
    pub fn build(time: f64) -> Self {
        let elapsed = time - T0;
        let step = Y0_MAX / N_STEPS as f64;
        let mut material_grid = Vec::with_capacity(N_STEPS + 1);
        let mut heights = Vec::with_capacity(N_STEPS + 1);

        let integrand = |height: f64| relative_growth_rate_y(height).powf(elapsed);
        let mut accumulated = 0.0;
        let mut previous = integrand(0.0);
        material_grid.push(0.0);
        heights.push(0.0);
        for index in 1..=N_STEPS {
            let material_height = index as f64 * step;
            let current = integrand(material_height);
            accumulated += 0.5 * (previous + current) * step;
            previous = current;
            material_grid.push(material_height);
            heights.push(accumulated);
        }
        VerticalMap {
            material_grid,
            heights,
        }
    }

    pub fn forward(&self, material_height: f64) -> f64 {
        interpolate(
            InterpolationAbscissae(&self.material_grid),
            InterpolationOrdinates(&self.heights),
            material_height,
        )
    }
    pub fn inverse(&self, height: f64) -> f64 {
        interpolate(
            InterpolationAbscissae(&self.heights),
            InterpolationOrdinates(&self.material_grid),
            height,
        )
    }
}

pub struct GrowthTransition<'a> {
    pub time_from: f64,
    pub time_to: f64,
    pub map_from: &'a VerticalMap,
    pub map_to: &'a VerticalMap,
}

pub fn propagate(point: Point, transition: &GrowthTransition) -> Point {
    let material_height = transition.map_from.inverse(point.y);
    let height = transition.map_to.forward(material_height);
    let x = point.x
        * relative_growth_rate_x(material_height).powf(transition.time_to - transition.time_from);
    Point::new(x, height)
}

pub fn unpropagate(point: Point, time: f64, map: &VerticalMap) -> Point {
    let material_height = map.inverse(point.y);
    let x = point.x / relative_growth_rate_x(material_height).powf(time - T0);
    Point::new(x, material_height)
}

pub fn reference_bounding_box(time: f64) -> BoundingBox {
    let scale_x = 1.30f64.powf(time);
    let scale_y = 1.25f64.powf(time);
    let half_width = 0.42 * scale_x * 1.15;
    BoundingBox {
        min_x: -half_width,
        max_x: half_width,
        min_y: 0.0,
        max_y: Y0_MAX * scale_y * 1.15,
    }
}
