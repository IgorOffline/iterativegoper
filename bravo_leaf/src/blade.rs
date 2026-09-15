use crate::constants::Y0_MAX;
use crate::geometry::Point;

pub fn leaf_half_width(height: f64) -> f64 {
    let taper = (std::f64::consts::PI * height).sin();
    0.42 * taper * (1.0 - 0.35 * height)
}

pub fn point_in_blade(point: Point) -> bool {
    (0.0..=Y0_MAX).contains(&point.y) && point.x.abs() <= leaf_half_width(point.y)
}

pub fn random_blade_point() -> Point {
    loop {
        let height = ::rand::random_range(0.0..Y0_MAX);
        let half_width = leaf_half_width(height);
        let candidate = Point::new(::rand::random_range(-half_width..=half_width), height);
        if point_in_blade(candidate) {
            return candidate;
        }
    }
}
