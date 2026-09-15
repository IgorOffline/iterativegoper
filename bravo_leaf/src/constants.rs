use macroquad::prelude::Color;

pub const SEED: u64 = 0x1EA_F_C0FFEE;

pub const T0: f64 = 0.0;
pub const T_MAX: f64 = 4.0;
pub const Y0_MAX: f64 = 1.0;
pub const N_STEPS: usize = 2048;

pub const BIRTH_DIST_SOURCE: f64 = 0.045;
pub const BIRTH_DIST_VEIN: f64 = 0.045;
pub const DART_ATTEMPTS: usize = 20_000;
pub const DART_ATTEMPTS_CYCLE: usize = 4_000;

pub const VEIN_GROWTH_STEP: f64 = 0.018;
pub const KILL_DISTANCE: f64 = 0.030;
pub const DEV_DT: f64 = 0.06;
pub const VENATION_SUBSTEPS: usize = 3;
pub const AUTO_STEP_INTERVAL: f64 = 0.05;

pub const MURRAY_EXPONENT: f64 = 3.0;
pub const TIP_RADIUS: f64 = 1.0;
pub const VEIN_WIDTH_PX: f32 = 1.4;

pub const MATERIAL_COLOR_BLUE: Color = Color::new(0.000, 0.478, 0.808, 1.0);
pub const BACKGROUND: Color = Color::new(0.08, 0.09, 0.11, 1.0);
