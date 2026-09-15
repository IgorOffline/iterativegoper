use macroquad::prelude::Color;

pub const SEEDS: [u64; 6] = [
    0x1EA_F_C0FFEE,
    0xB105_EED5_0DDBA11,
    0x0FF1CE_C0DE_CAFE,
    0xACE1_DEAD_BEEF_F00D,
    0xDEC0DED_1EAF_D00D,
    0xFACADE_D15EA5ED,
];

pub const MATERIAL_PALETTE_V3_SHADE_500: [Color; 6] = [
    Color::new(0.957, 0.263, 0.212, 1.0),
    Color::new(0.298, 0.686, 0.314, 1.0),
    Color::new(0.129, 0.588, 0.953, 1.0),
    Color::new(0.612, 0.153, 0.690, 1.0),
    Color::new(0.376, 0.490, 0.545, 1.0),
    Color::new(1.000, 0.596, 0.000, 1.0),
];

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

pub const SUBJECT_LIFESPAN: f64 = 6.0;
pub const SUBJECT_SPAWN_INTERVAL: f64 = 4.0;
pub const SUBJECT_FADE: f64 = 1.0;

pub const BACKGROUND: Color = Color::new(0.08, 0.09, 0.11, 1.0);
