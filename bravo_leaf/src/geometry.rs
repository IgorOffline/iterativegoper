#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub const fn new(x: f64, y: f64) -> Self {
        Point { x, y }
    }
}

pub fn distance_squared(first: Point, second: Point) -> f64 {
    let delta_x = first.x - second.x;
    let delta_y = first.y - second.y;
    delta_x * delta_x + delta_y * delta_y
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub min_x: f64,
    pub max_x: f64,
    pub min_y: f64,
    pub max_y: f64,
}
