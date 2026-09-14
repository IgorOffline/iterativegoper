#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PointX {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PointY {
    pub ada: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Point {
    x: PointX,
    y: PointY,
}

impl Point {
    fn new(x: PointX, y: PointY) -> Self {
        Point { x, y }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GrowthModelLeft {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GrowthModelRight {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct VerticalMapY0 {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct VerticalMapValue {
    pub ada: f64,
}

struct VerticalMap {
    y0_grid: Vec<VerticalMapY0>,
    values: Vec<VerticalMapValue>,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct VerticalMapTau {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct VerticalMapY0Max {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct VerticalMapNSteps {
    pub ada: usize,
}

impl VerticalMap {
    fn build<RY>(
        rerg_y: RY,
        tau: VerticalMapTau,
        y0_max: VerticalMapY0Max,
        n_steps: VerticalMapNSteps,
    ) -> Self
    where
        RY: Fn(GrowthModelRight) -> GrowthModelRight,
    {
        let tau = tau.ada;
        let y0_max = y0_max.ada;
        let n_steps = n_steps.ada;

        assert!(y0_max > 0.0, "y0_max must be positive");
        assert!(n_steps >= 1, "need at least one integration step");

        let h = y0_max / n_steps as f64;
        let mut y0_grid = Vec::with_capacity(n_steps + 1);
        let mut values = Vec::with_capacity(n_steps + 1);

        let f = |y: f64| rerg_y(GrowthModelRight { ada: y }).ada.powf(tau);

        let mut acc = 0.0;
        let mut prev_f = f(0.0);
        y0_grid.push(VerticalMapY0 { ada: 0.0 });
        values.push(VerticalMapValue { ada: 0.0 });

        for i in 1..=n_steps {
            let y = i as f64 * h;
            let cur_f = f(y);
            acc += 0.5 * (prev_f + cur_f) * h;
            prev_f = cur_f;
            y0_grid.push(VerticalMapY0 { ada: y });
            values.push(VerticalMapValue { ada: acc });
        }

        VerticalMap { y0_grid, values }
    }

    fn forward(&self, y0: f64) -> f64 {
        interp(
            &self.y0_grid,
            &self.values,
            InterpX { ada: y0 },
            |g| InterpLeft { ada: g.ada },
            |v| InterpRight { ada: v.ada },
        )
        .ada
    }

    fn inverse(&self, y: f64) -> f64 {
        interp(
            &self.values,
            &self.y0_grid,
            InterpX { ada: y },
            |v| InterpLeft { ada: v.ada },
            |g| InterpRight { ada: g.ada },
        )
        .ada
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct InterpLeft {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct InterpRight {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct InterpX {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct InterpRetVal {
    pub ada: f64,
}

fn interp<L, R, IL, IR>(xs: &[L], ys: &[R], x: InterpX, xa: IL, ya: IR) -> InterpRetVal
where
    IL: Fn(&L) -> InterpLeft,
    IR: Fn(&R) -> InterpRight,
{
    debug_assert_eq!(xs.len(), ys.len());
    let x = x.ada;
    let n = xs.len();
    if x <= xa(&xs[0]).ada {
        return InterpRetVal {
            ada: ya(&ys[0]).ada,
        };
    }
    if x >= xa(&xs[n - 1]).ada {
        return InterpRetVal {
            ada: ya(&ys[n - 1]).ada,
        };
    }
    let mut lo = 0usize;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xa(&xs[mid]).ada <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (x - xa(&xs[lo]).ada) / (xa(&xs[hi]).ada - xa(&xs[lo]).ada);
    InterpRetVal {
        ada: ya(&ys[lo]).ada + t * (ya(&ys[hi]).ada - ya(&ys[lo]).ada),
    }
}

struct GrowthModel<FL, FR>
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    rerg_x: FL,
    rerg_y: FR,
    t0: f64,
}

impl<FL, FR> GrowthModel<FL, FR>
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    fn new(rerg_x: FL, rerg_y: FR, t0: f64) -> Self {
        GrowthModel { rerg_x, rerg_y, t0 }
    }

    fn rerg_x_at(&self, y0: f64) -> f64 {
        (self.rerg_x)(GrowthModelLeft { ada: y0 }).ada
    }

    fn vertical_map(&self, t: f64, y0_max: f64, n_steps: usize) -> VerticalMap {
        VerticalMap::build(
            &self.rerg_y,
            VerticalMapTau { ada: t - self.t0 },
            VerticalMapY0Max { ada: y0_max },
            VerticalMapNSteps { ada: n_steps },
        )
    }

    fn propagate_point(
        &self,
        p: Point,
        t1: f64,
        t2: f64,
        map_t1: &VerticalMap,
        map_t2: &VerticalMap,
    ) -> Point {
        let y0 = map_t1.inverse(p.y.ada);
        let y2 = map_t2.forward(y0);
        let x2 = p.x.ada * self.rerg_x_at(y0).powf(t2 - t1);
        Point::new(PointX { ada: x2 }, PointY { ada: y2 })
    }

    #[allow(dead_code)]
    fn propagate_all(
        &self,
        points: &[Point],
        t1: f64,
        t2: f64,
        y0_max: f64,
        n_steps: usize,
    ) -> Vec<Point> {
        let map_t1 = self.vertical_map(t1, y0_max, n_steps);
        let map_t2 = self.vertical_map(t2, y0_max, n_steps);
        points
            .iter()
            .map(|&p| self.propagate_point(p, t1, t2, &map_t1, &map_t2))
            .collect()
    }
}

use macroquad::prelude::*;

const T0: f64 = 0.0;
const T_MAX: f64 = 4.0;
const Y0_MAX: f64 = 1.0;
const N_STEPS: usize = 2048;

const BIRTH_DIST_SOURCE: f64 = 0.045;
const BIRTH_DIST_VEIN: f64 = 0.045;
const DART_ATTEMPTS: usize = 20000;

const VEIN_GROWTH_STEP: f64 = 0.018;
const KILL_DISTANCE: f64 = 0.030;
const AUTO_STEP_INTERVAL: f64 = 0.05;

fn leaf_half_width(y: f64) -> f64 {
    let s = (std::f64::consts::PI * y).sin();
    0.42 * s * (1.0 - 0.35 * y)
}

fn point_in_blade(x: f64, y: f64) -> bool {
    (0.0..=Y0_MAX).contains(&y) && x.abs() <= leaf_half_width(y)
}

fn reference_leaf() -> Vec<Point> {
    let n = 120usize;
    let mut outline = Vec::with_capacity(2 * n + 1);
    for i in 0..=n {
        let y = Y0_MAX * i as f64 / n as f64;
        outline.push(Point::new(
            PointX {
                ada: leaf_half_width(y),
            },
            PointY { ada: y },
        ));
    }
    for i in (0..=n).rev() {
        let y = Y0_MAX * i as f64 / n as f64;
        outline.push(Point::new(
            PointX {
                ada: -leaf_half_width(y),
            },
            PointY { ada: y },
        ));
    }

    outline
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct VeinNodeId {
    pub ada: usize,
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct VeinNode {
    pos: Point,
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct VeinEdge {
    from: VeinNodeId,
    to: VeinNodeId,
}

#[derive(Debug, Clone, PartialEq)]
struct VeinGraph {
    nodes: Vec<VeinNode>,
    edges: Vec<VeinEdge>,
}

impl VeinGraph {
    fn new() -> Self {
        VeinGraph {
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    fn add_node(&mut self, pos: Point) -> VeinNodeId {
        let id = VeinNodeId {
            ada: self.nodes.len(),
        };
        self.nodes.push(VeinNode { pos });
        id
    }

    fn add_edge(&mut self, from: VeinNodeId, to: VeinNodeId) {
        self.edges.push(VeinEdge { from, to });
    }
}

fn seed_graph() -> VeinGraph {
    let mut g = VeinGraph::new();
    g.add_node(Point::new(PointX { ada: 0.0 }, PointY { ada: 0.0 }));
    g
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct Source {
    pos: Point,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct BirthDistanceSource {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct BirthDistanceVein {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct DartAttempts {
    pub ada: usize,
}

fn dist2(a: &Point, b: &Point) -> f64 {
    let dx = a.x.ada - b.x.ada;
    let dy = a.y.ada - b.y.ada;
    dx * dx + dy * dy
}

fn throw_darts(
    graph: &VeinGraph,
    b_s: BirthDistanceSource,
    b_v: BirthDistanceVein,
    attempts: DartAttempts,
) -> Vec<Source> {
    let b_s = b_s.ada;
    let b_v = b_v.ada;
    let b_s2 = b_s * b_s;
    let b_v2 = b_v * b_v;

    let mut sources: Vec<Source> = Vec::new();
    for _ in 0..attempts.ada {
        let y = rand::gen_range(0.0f64, Y0_MAX);
        let hw = leaf_half_width(y);
        let x = rand::gen_range(-hw, hw);
        if !point_in_blade(x, y) {
            continue;
        }
        let cand = Point::new(PointX { ada: x }, PointY { ada: y });

        let far_from_sources = sources.iter().all(|s| dist2(&s.pos, &cand) >= b_s2);
        if !far_from_sources {
            continue;
        }
        let far_from_veins = graph.nodes.iter().all(|nd| dist2(&nd.pos, &cand) >= b_v2);
        if !far_from_veins {
            continue;
        }

        sources.push(Source { pos: cand });
    }
    sources
}

#[allow(dead_code)]
fn propagate_sources<FL, FR>(
    model: &GrowthModel<FL, FR>,
    sources: &[Source],
    map_t0: &VerticalMap,
    map_t: &VerticalMap,
    t: f64,
) -> Vec<Source>
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    sources
        .iter()
        .map(|s| Source {
            pos: model.propagate_point(s.pos, T0, t, map_t0, map_t),
        })
        .collect()
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GrowthStep {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct KillDistance {
    pub ada: f64,
}

fn nearest_node(graph: &VeinGraph, p: &Point) -> Option<VeinNodeId> {
    let mut best: Option<(usize, f64)> = None;
    for (i, nd) in graph.nodes.iter().enumerate() {
        let d = dist2(&nd.pos, p);
        match best {
            Some((_, bd)) if d >= bd => {}
            _ => best = Some((i, d)),
        }
    }
    best.map(|(i, _)| VeinNodeId { ada: i })
}

fn venation_step(
    graph: &mut VeinGraph,
    sources: &mut Vec<Source>,
    d: GrowthStep,
    dk: KillDistance,
) -> usize {
    let d = d.ada;
    let dk2 = dk.ada * dk.ada;

    if sources.is_empty() || graph.nodes.is_empty() {
        return 0;
    }

    let mut influence: Vec<(f64, f64)> = vec![(0.0, 0.0); graph.nodes.len()];
    let mut has_influence = vec![false; graph.nodes.len()];

    for s in sources.iter() {
        if let Some(v) = nearest_node(graph, &s.pos) {
            let vpos = graph.nodes[v.ada].pos;
            let dx = s.pos.x.ada - vpos.x.ada;
            let dy = s.pos.y.ada - vpos.y.ada;
            let len = (dx * dx + dy * dy).sqrt();
            if len > 1e-12 {
                influence[v.ada].0 += dx / len;
                influence[v.ada].1 += dy / len;
                has_influence[v.ada] = true;
            }
        }
    }

    let mut spawned = 0usize;
    let original_len = graph.nodes.len();
    for v in 0..original_len {
        if !has_influence[v] {
            continue;
        }
        let (ax, ay) = influence[v];
        let nlen = (ax * ax + ay * ay).sqrt();
        if nlen <= 1e-12 {
            continue;
        }
        let vpos = graph.nodes[v].pos;
        let nx = ax / nlen;
        let ny = ay / nlen;
        let new_pos = Point::new(
            PointX {
                ada: vpos.x.ada + d * nx,
            },
            PointY {
                ada: vpos.y.ada + d * ny,
            },
        );
        let new_id = graph.add_node(new_pos);
        graph.add_edge(VeinNodeId { ada: v }, new_id);
        spawned += 1;
    }

    sources.retain(|s| graph.nodes.iter().all(|nd| dist2(&nd.pos, &s.pos) > dk2));

    spawned
}

#[allow(dead_code)]
fn propagate_graph<FL, FR>(
    model: &GrowthModel<FL, FR>,
    graph: &VeinGraph,
    map_t0: &VerticalMap,
    map_t: &VerticalMap,
    t: f64,
) -> VeinGraph
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    let nodes = graph
        .nodes
        .iter()
        .map(|nd| VeinNode {
            pos: model.propagate_point(nd.pos, T0, t, map_t0, map_t),
        })
        .collect();
    VeinGraph {
        nodes,
        edges: graph.edges.clone(),
    }
}

fn propagate_from_ref<FL, FR>(
    model: &GrowthModel<FL, FR>,
    pts: &[Point],
    map_t0: &VerticalMap,
    map_t: &VerticalMap,
    t: f64,
) -> Vec<Point>
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    pts.iter()
        .map(|&p| model.propagate_point(p, T0, t, map_t0, map_t))
        .collect()
}

struct Fit {
    scale: f32,
    ox: f32,
    oy: f32,
}

fn fit_world(bbox: (f64, f64, f64, f64)) -> Fit {
    let (minx, maxx, miny, maxy) = bbox;
    let margin = 60.0f32;
    let w = screen_width() - 2.0 * margin;
    let h = screen_height() - 2.0 * margin - 40.0;
    let span_x = (maxx - minx).max(1e-6) as f32;
    let span_y = (maxy - miny).max(1e-6) as f32;
    let scale = (w / span_x).min(h / span_y);
    let cx = 0.5 * (minx + maxx) as f32;
    let ox = screen_width() * 0.5 - cx * scale;
    let oy = screen_height() - margin;
    Fit { scale, ox, oy }
}

fn to_screen(p: &Point, fit: &Fit) -> (f32, f32) {
    let sx = fit.ox + p.x.ada as f32 * fit.scale;
    let sy = fit.oy - p.y.ada as f32 * fit.scale;
    (sx, sy)
}

fn draw_polyline(pts: &[Point], fit: &Fit, color: Color, thick: f32, closed: bool) {
    if pts.len() < 2 {
        return;
    }
    for w in pts.windows(2) {
        let (x1, y1) = to_screen(&w[0], fit);
        let (x2, y2) = to_screen(&w[1], fit);
        draw_line(x1, y1, x2, y2, thick, color);
    }
    if closed {
        let (x1, y1) = to_screen(&pts[pts.len() - 1], fit);
        let (x2, y2) = to_screen(&pts[0], fit);
        draw_line(x1, y1, x2, y2, thick, color);
    }
}

fn draw_vein_graph(graph: &VeinGraph, fit: &Fit, edge_color: Color, node_color: Color) {
    for e in &graph.edges {
        let a = &graph.nodes[e.from.ada].pos;
        let b = &graph.nodes[e.to.ada].pos;
        let (x1, y1) = to_screen(a, fit);
        let (x2, y2) = to_screen(b, fit);
        draw_line(x1, y1, x2, y2, 2.0, edge_color);
    }
    for (i, nd) in graph.nodes.iter().enumerate() {
        let (sx, sy) = to_screen(&nd.pos, fit);
        if i == 0 {
            draw_circle_lines(sx, sy, 9.0, 2.0, node_color);
            draw_circle(sx, sy, 4.5, node_color);
        } else {
            draw_circle(sx, sy, 3.5, node_color);
        }
    }
}

fn draw_sources(sources: &[Source], fit: &Fit, color: Color) {
    for s in sources {
        let (sx, sy) = to_screen(&s.pos, fit);
        draw_circle(sx, sy, 2.5, color);
    }
}

fn ref_bbox(t: f64) -> (f64, f64, f64, f64) {
    let sx = 1.30f64.powf(t);
    let sy = 1.25f64.powf(t);
    let hx = 0.42 * sx * 1.15;
    (-hx, hx, 0.0, Y0_MAX * sy * 1.15)
}

#[macroquad::main("Nonuniform leaf growth")]
async fn main() {
    let rerg_x = |l: GrowthModelLeft| GrowthModelLeft {
        ada: 1.30 - 0.20 * l.ada,
    };
    let rerg_y = |r: GrowthModelRight| GrowthModelRight {
        ada: 1.25 - 0.15 * r.ada,
    };
    let model = GrowthModel::new(rerg_x, rerg_y, T0);

    let outline_ref = reference_leaf();
    let map_t0 = model.vertical_map(T0, Y0_MAX, N_STEPS);

    let dart = |graph: &VeinGraph| {
        throw_darts(
            graph,
            BirthDistanceSource {
                ada: BIRTH_DIST_SOURCE,
            },
            BirthDistanceVein {
                ada: BIRTH_DIST_VEIN,
            },
            DartAttempts { ada: DART_ATTEMPTS },
        )
    };

    // Step 3 is venation-first: time is frozen at T0 and veins grow in
    // reference space. graph and sources are persistent, iterating state.
    let mut graph = seed_graph();
    let mut sources = dart(&graph);
    let initial_sources = sources.len();
    let mut iteration = 0usize;

    let t = T0;
    let map_t = model.vertical_map(t, Y0_MAX, N_STEPS);
    let outline = propagate_from_ref(&model, &outline_ref, &map_t0, &map_t, t);

    let mut auto = false;
    let mut auto_accum = 0.0f64;
    let mut stall = 0usize;
    let stall_limit = 400usize;

    let leaf_green = Color::new(0.20, 0.55, 0.25, 1.0);
    let vein_green = Color::new(0.35, 0.70, 0.40, 0.9);
    let node_yellow = Color::new(0.95, 0.85, 0.30, 1.0);
    let source_blue = Color::new(0.45, 0.70, 0.95, 0.9);

    let do_step = |graph: &mut VeinGraph, sources: &mut Vec<Source>| {
        venation_step(
            graph,
            sources,
            GrowthStep {
                ada: VEIN_GROWTH_STEP,
            },
            KillDistance { ada: KILL_DISTANCE },
        )
    };

    loop {
        if is_key_pressed(KeyCode::Space) {
            auto = !auto;
        }
        if is_key_pressed(KeyCode::N) {
            do_step(&mut graph, &mut sources);
            iteration += 1;
        }
        if is_key_pressed(KeyCode::R) {
            graph = seed_graph();
            sources = dart(&graph);
            iteration = 0;
            auto = false;
            stall = 0;
        }
        if is_key_pressed(KeyCode::S) {
            sources = dart(&graph);
            stall = 0;
        }
        if auto && !sources.is_empty() && stall < stall_limit {
            auto_accum += get_frame_time() as f64;
            while auto_accum >= AUTO_STEP_INTERVAL {
                auto_accum -= AUTO_STEP_INTERVAL;
                let before = sources.len();
                do_step(&mut graph, &mut sources);
                iteration += 1;
                if sources.len() < before {
                    stall = 0;
                } else {
                    stall += 1;
                }
            }
            if stall >= stall_limit {
                auto = false;
            }
        }

        let fit = fit_world(ref_bbox(T_MAX));

        clear_background(Color::new(0.08, 0.09, 0.11, 1.0));

        draw_polyline(&outline, &fit, leaf_green, 2.5, true);

        for &p in outline.iter().step_by(20) {
            let (sx, sy) = to_screen(&p, &fit);
            draw_circle(sx, sy, 2.5, leaf_green);
        }

        draw_sources(&sources, &fit, source_blue);
        draw_vein_graph(&graph, &fit, vein_green, node_yellow);

        let hud_title = format!(
            "Step 3: open venation (Eq. 1)   iteration = {iteration}   [auto: {}]",
            if auto { "on" } else { "off" }
        );
        draw_text(&hud_title, 20.0, 30.0, 26.0, WHITE);
        let hud_veins = format!(
            "vein nodes = {}  edges = {}   |   D = {VEIN_GROWTH_STEP:.3}  d_k = {KILL_DISTANCE:.3}",
            graph.nodes.len(),
            graph.edges.len()
        );
        draw_text(&hud_veins, 20.0, 56.0, 22.0, node_yellow);
        let status = if sources.is_empty() {
            "  [complete: all sources consumed]"
        } else if stall >= stall_limit {
            "  [settled: remaining sources orphaned -- open pattern]"
        } else {
            ""
        };
        let hud_sources = format!(
            "auxin sources: {} remaining of {} (consumed = {}){status}",
            sources.len(),
            initial_sources,
            initial_sources.saturating_sub(sources.len())
        );
        draw_text(&hud_sources, 20.0, 82.0, 22.0, source_blue);
        draw_text(
            "[n] one step   [space] auto   [r] reset   [s] re-throw sources   (t frozen at t0 -- growth returns in Step 4)",
            20.0,
            screen_height() - 16.0,
            20.0,
            GRAY,
        );

        next_frame().await
    }
}
