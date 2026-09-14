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

    fn unpropagate_point(&self, p: Point, t: f64, map_t: &VerticalMap) -> Point {
        let y0 = map_t.inverse(p.y.ada);
        let x0 = p.x.ada / self.rerg_x_at(y0).powf(t - self.t0);
        Point::new(PointX { ada: x0 }, PointY { ada: y0 })
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

const DEV_DT: f64 = 0.06;
const VENATION_SUBSTEPS: usize = 3;
const DART_ATTEMPTS_CYCLE: usize = 4000;

const MURRAY_EXPONENT: f64 = 3.0;
const TIP_RADIUS: f64 = 1.0;
const VEIN_WIDTH_PX: f32 = 1.4;

const MERGE_DISTANCE: f64 = 0.022;

const SPATIAL_INDEX_MIN_NODES: usize = 1500;

const TOOTH_COUNT: f64 = 11.0;
const TOOTH_DEPTH: f64 = 0.06;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum LeafForm {
    SmoothPinnate,
    ToothedPinnate,
    ToothedActinodromous,
}

static LEAF_FORM: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);

fn set_leaf_form(f: LeafForm) {
    let v = match f {
        LeafForm::SmoothPinnate => 0,
        LeafForm::ToothedPinnate => 1,
        LeafForm::ToothedActinodromous => 2,
    };
    LEAF_FORM.store(v, std::sync::atomic::Ordering::Relaxed);
}

fn leaf_form() -> LeafForm {
    match LEAF_FORM.load(std::sync::atomic::Ordering::Relaxed) {
        1 => LeafForm::ToothedPinnate,
        2 => LeafForm::ToothedActinodromous,
        _ => LeafForm::SmoothPinnate,
    }
}

fn triangular_wave(u: f64, period: f64) -> f64 {
    let phase = (u / period).rem_euclid(1.0);
    1.0 - (2.0 * phase - 1.0).abs()
}

fn tooth_offset(y: f64) -> f64 {
    let taper = (std::f64::consts::PI * y).sin();
    let primary = triangular_wave(y, 1.0 / TOOTH_COUNT);
    let secondary = 0.4 * triangular_wave(y, 1.0 / (2.0 * TOOTH_COUNT));
    TOOTH_DEPTH * taper * (primary + secondary)
}

fn leaf_half_width(y: f64) -> f64 {
    let s = (std::f64::consts::PI * y).sin();
    let base = 0.42 * s * (1.0 - 0.35 * y);
    match leaf_form() {
        LeafForm::SmoothPinnate => base,
        LeafForm::ToothedPinnate | LeafForm::ToothedActinodromous => base + tooth_offset(y),
    }
}

fn point_in_blade(x: f64, y: f64) -> bool {
    (0.0..=Y0_MAX).contains(&y) && x.abs() <= leaf_half_width(y)
}

fn reference_leaf() -> Vec<Point> {
    let n = 400usize;
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
    let root = g.add_node(Point::new(PointX { ada: 0.0 }, PointY { ada: 0.0 }));

    if leaf_form() == LeafForm::ToothedActinodromous {
        let primaries = 5usize;
        let stub = 0.05f64;
        for k in 0..primaries {
            let frac = (k as f64 + 0.5) / primaries as f64;
            let angle = (frac - 0.5) * std::f64::consts::PI * 0.9;
            let px = stub * angle.sin();
            let py = (stub * angle.cos()).max(1e-3);
            let id = g.add_node(Point::new(PointX { ada: px }, PointY { ada: py }));
            g.add_edge(root, id);
        }
    }
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

use spade::{DelaunayTriangulation, HasPosition, Point2, Triangulation};

struct IndexedVertex {
    point: Point2<f64>,
    node: usize,
}

impl HasPosition for IndexedVertex {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        self.point
    }
}

struct NodeIndex {
    tri: DelaunayTriangulation<IndexedVertex>,
    ok: bool,
}

impl NodeIndex {
    fn build(graph: &VeinGraph) -> Self {
        let mut tri: DelaunayTriangulation<IndexedVertex> = DelaunayTriangulation::new();
        let mut ok = true;
        for (i, nd) in graph.nodes.iter().enumerate() {
            let v = IndexedVertex {
                point: Point2::new(nd.pos.x.ada, nd.pos.y.ada),
                node: i,
            };
            if tri.insert(v).is_err() {
                ok = false;
                break;
            }
        }
        if tri.num_vertices() < 2 {
            ok = false;
        }
        NodeIndex { tri, ok }
    }

    fn nearest(&self, graph: &VeinGraph, p: &Point) -> Option<usize> {
        if self.ok
            && let Some(vh) = self.tri.nearest_neighbor(Point2::new(p.x.ada, p.y.ada))
        {
            return Some(vh.data().node);
        }
        nearest_node_linear(graph, p)
    }
}

fn nearest_node_linear(graph: &VeinGraph, p: &Point) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (i, nd) in graph.nodes.iter().enumerate() {
        let d = dist2(&nd.pos, p);
        match best {
            Some((_, bd)) if d >= bd => {}
            _ => best = Some((i, d)),
        }
    }
    best.map(|(i, _)| i)
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

fn point_in_current_blade<FL, FR>(
    model: &GrowthModel<FL, FR>,
    p: Point,
    t: f64,
    map_t: &VerticalMap,
) -> bool
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    let r = model.unpropagate_point(p, t, map_t);
    point_in_blade(r.x.ada, r.y.ada)
}

#[allow(clippy::too_many_arguments)]
fn throw_darts_current<FL, FR>(
    model: &GrowthModel<FL, FR>,
    t: f64,
    map_t: &VerticalMap,
    graph: &VeinGraph,
    existing: &[Source],
    bbox: (f64, f64, f64, f64),
    b_s: BirthDistanceSource,
    b_v: BirthDistanceVein,
    attempts: DartAttempts,
) -> Vec<Source>
where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    let b_s2 = b_s.ada * b_s.ada;
    let b_v2 = b_v.ada * b_v.ada;
    let (minx, maxx, miny, maxy) = bbox;

    let mut fresh: Vec<Source> = Vec::new();
    for _ in 0..attempts.ada {
        let x = rand::gen_range(minx, maxx);
        let y = rand::gen_range(miny, maxy);
        let cand = Point::new(PointX { ada: x }, PointY { ada: y });
        if !point_in_current_blade(model, cand, t, map_t) {
            continue;
        }
        if existing.iter().any(|s| dist2(&s.pos, &cand) < b_s2) {
            continue;
        }
        if fresh.iter().any(|s| dist2(&s.pos, &cand) < b_s2) {
            continue;
        }
        if graph.nodes.iter().any(|nd| dist2(&nd.pos, &cand) < b_v2) {
            continue;
        }
        fresh.push(Source { pos: cand });
    }
    fresh
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GrowthStep {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct KillDistance {
    pub ada: f64,
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

    let index = if graph.nodes.len() >= SPATIAL_INDEX_MIN_NODES {
        Some(NodeIndex::build(graph))
    } else {
        None
    };
    let nearest = |g: &VeinGraph, p: &Point| -> Option<usize> {
        match &index {
            Some(idx) => idx.nearest(g, p),
            None => nearest_node_linear(g, p),
        }
    };

    let mut influence: Vec<(f64, f64)> = vec![(0.0, 0.0); graph.nodes.len()];
    let mut has_influence = vec![false; graph.nodes.len()];

    for s in sources.iter() {
        if let Some(v) = nearest(graph, &s.pos) {
            let vpos = graph.nodes[v].pos;
            let dx = s.pos.x.ada - vpos.x.ada;
            let dy = s.pos.y.ada - vpos.y.ada;
            let len = (dx * dx + dy * dy).sqrt();
            if len > 1e-12 {
                influence[v].0 += dx / len;
                influence[v].1 += dy / len;
                has_influence[v] = true;
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

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum VenationMode {
    Open,
    Closed,
}

fn source_influences_node(graph: &VeinGraph, index: &NodeIndex, s: &Point, v: usize) -> bool {
    let vpos = graph.nodes[v].pos;
    let d_sv = dist2(s, &vpos);

    if let Some(w0) = index.nearest(graph, s)
        && w0 != v
    {
        let nd = graph.nodes[w0].pos;
        if dist2(s, &nd) < d_sv && dist2(&vpos, &nd) < d_sv {
            return false;
        }
    }

    for (w, nd) in graph.nodes.iter().enumerate() {
        if w == v {
            continue;
        }
        let d_sw = dist2(s, &nd.pos);
        let d_vw = dist2(&vpos, &nd.pos);
        if d_sw < d_sv && d_vw < d_sv {
            return false;
        }
    }
    true
}

fn venation_step_closed(
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

    let index = NodeIndex::build(graph);

    let mut influence: Vec<(f64, f64)> = vec![(0.0, 0.0); graph.nodes.len()];
    let mut has_influence = vec![false; graph.nodes.len()];

    for s in sources.iter() {
        for v in 0..graph.nodes.len() {
            if !source_influences_node(graph, &index, &s.pos, v) {
                continue;
            }
            let vpos = graph.nodes[v].pos;
            let dx = s.pos.x.ada - vpos.x.ada;
            let dy = s.pos.y.ada - vpos.y.ada;
            let len = (dx * dx + dy * dy).sqrt();
            if len > 1e-12 {
                influence[v].0 += dx / len;
                influence[v].1 += dy / len;
                has_influence[v] = true;
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

    close_loops(graph, original_len);

    sources.retain(|s| graph.nodes.iter().all(|nd| dist2(&nd.pos, &s.pos) > dk2));

    spawned
}

fn close_loops(graph: &mut VeinGraph, first_new: usize) {
    let merge2 = MERGE_DISTANCE * MERGE_DISTANCE;
    let mut new_edges: Vec<(usize, usize)> = Vec::new();
    for nw in first_new..graph.nodes.len() {
        let np = graph.nodes[nw].pos;
        let parent = graph
            .edges
            .iter()
            .find(|e| e.to.ada == nw)
            .map(|e| e.from.ada);
        let mut best: Option<(usize, f64)> = None;
        for old in 0..first_new {
            if Some(old) == parent {
                continue;
            }
            let d = dist2(&np, &graph.nodes[old].pos);
            if d <= merge2 {
                match best {
                    Some((_, bd)) if d >= bd => {}
                    _ => best = Some((old, d)),
                }
            }
        }
        if let Some((old, _)) = best {
            new_edges.push((nw, old));
        }
    }
    for (nw, old) in new_edges {
        graph.add_edge(VeinNodeId { ada: nw }, VeinNodeId { ada: old });
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct MurrayExponent {
    pub ada: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct TipRadius {
    pub ada: f64,
}

fn murray_radii(graph: &VeinGraph, n: MurrayExponent, r0: TipRadius) -> Vec<f64> {
    let n = n.ada;
    let r0 = r0.ada;
    let count = graph.nodes.len();

    let mut child_pow_sum = vec![0.0f64; count];
    let mut has_child = vec![false; count];
    let mut parent = vec![None::<usize>; count];
    for e in &graph.edges {
        if e.from.ada < e.to.ada {
            parent[e.to.ada] = Some(e.from.ada);
        }
    }

    let mut radius = vec![r0; count];
    for v in (0..count).rev() {
        radius[v] = if has_child[v] {
            child_pow_sum[v].powf(1.0 / n)
        } else {
            r0
        };
        if let Some(p) = parent[v] {
            child_pow_sum[p] += radius[v].powf(n);
            has_child[p] = true;
        }
    }

    radius
}

fn displace_graph<FL, FR>(
    model: &GrowthModel<FL, FR>,
    graph: &mut VeinGraph,
    t1: f64,
    t2: f64,
    map_t1: &VerticalMap,
    map_t2: &VerticalMap,
) where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    for nd in graph.nodes.iter_mut() {
        nd.pos = model.propagate_point(nd.pos, t1, t2, map_t1, map_t2);
    }
}

fn displace_sources<FL, FR>(
    model: &GrowthModel<FL, FR>,
    sources: &mut [Source],
    t1: f64,
    t2: f64,
    map_t1: &VerticalMap,
    map_t2: &VerticalMap,
) where
    FL: Fn(GrowthModelLeft) -> GrowthModelLeft,
    FR: Fn(GrowthModelRight) -> GrowthModelRight,
{
    for s in sources.iter_mut() {
        s.pos = model.propagate_point(s.pos, t1, t2, map_t1, map_t2);
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

fn draw_vein_graph(
    graph: &VeinGraph,
    radii: &[f64],
    fit: &Fit,
    edge_color: Color,
    node_color: Color,
) {
    for e in &graph.edges {
        let a = &graph.nodes[e.from.ada].pos;
        let b = &graph.nodes[e.to.ada].pos;
        let (x1, y1) = to_screen(a, fit);
        let (x2, y2) = to_screen(b, fit);
        let w = if e.from.ada < e.to.ada {
            (radii[e.to.ada] as f32 * VEIN_WIDTH_PX).max(1.0)
        } else {
            1.0
        };
        draw_line(x1, y1, x2, y2, w, edge_color);
    }
    if let Some(nd) = graph.nodes.first() {
        let (sx, sy) = to_screen(&nd.pos, fit);
        draw_circle_lines(sx, sy, 9.0, 2.0, node_color);
        draw_circle(sx, sy, 4.5, node_color);
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

    set_leaf_form(LeafForm::SmoothPinnate);
    let mut form = LeafForm::SmoothPinnate;

    let mut outline_ref = reference_leaf();
    let map_t0 = model.vertical_map(T0, Y0_MAX, N_STEPS);

    let mut graph = seed_graph();
    let mut sources: Vec<Source> = throw_darts(
        &graph,
        BirthDistanceSource {
            ada: BIRTH_DIST_SOURCE,
        },
        BirthDistanceVein {
            ada: BIRTH_DIST_VEIN,
        },
        DartAttempts { ada: DART_ATTEMPTS },
    );
    let mut t = T0;
    let mut cycle = 0usize;
    let mut mode = VenationMode::Open;

    let mut auto = true;
    let mut auto_accum = 0.0f64;

    let leaf_green = Color::new(0.20, 0.55, 0.25, 1.0);
    let vein_green = Color::new(0.35, 0.70, 0.40, 0.9);
    let node_yellow = Color::new(0.95, 0.85, 0.30, 1.0);
    let source_blue = Color::new(0.45, 0.70, 0.95, 0.9);

    let grow_cycle =
        |graph: &mut VeinGraph, sources: &mut Vec<Source>, t: &mut f64, mode: VenationMode| {
            if *t >= T_MAX {
                return;
            }
            let t1 = *t;
            let t2 = (t1 + DEV_DT).min(T_MAX);
            let map_t1 = model.vertical_map(t1, Y0_MAX, N_STEPS);
            let map_t2 = model.vertical_map(t2, Y0_MAX, N_STEPS);

            displace_graph(&model, graph, t1, t2, &map_t1, &map_t2);
            displace_sources(&model, sources, t1, t2, &map_t1, &map_t2);
            *t = t2;

            let fresh = throw_darts_current(
                &model,
                t2,
                &map_t2,
                graph,
                sources,
                ref_bbox(t2),
                BirthDistanceSource {
                    ada: BIRTH_DIST_SOURCE,
                },
                BirthDistanceVein {
                    ada: BIRTH_DIST_VEIN,
                },
                DartAttempts {
                    ada: DART_ATTEMPTS_CYCLE,
                },
            );
            sources.extend(fresh);

            for _ in 0..VENATION_SUBSTEPS {
                match mode {
                    VenationMode::Open => venation_step(
                        graph,
                        sources,
                        GrowthStep {
                            ada: VEIN_GROWTH_STEP,
                        },
                        KillDistance { ada: KILL_DISTANCE },
                    ),
                    VenationMode::Closed => venation_step_closed(
                        graph,
                        sources,
                        GrowthStep {
                            ada: VEIN_GROWTH_STEP,
                        },
                        KillDistance { ada: KILL_DISTANCE },
                    ),
                };
            }
        };

    loop {
        if is_key_pressed(KeyCode::Space) {
            auto = !auto;
        }
        if is_key_pressed(KeyCode::C) {
            grow_cycle(&mut graph, &mut sources, &mut t, mode);
            cycle += 1;
        }
        if is_key_pressed(KeyCode::O) {
            mode = match mode {
                VenationMode::Open => VenationMode::Closed,
                VenationMode::Closed => VenationMode::Open,
            };
            graph = seed_graph();
            sources = throw_darts(
                &graph,
                BirthDistanceSource {
                    ada: BIRTH_DIST_SOURCE,
                },
                BirthDistanceVein {
                    ada: BIRTH_DIST_VEIN,
                },
                DartAttempts { ada: DART_ATTEMPTS },
            );
            t = T0;
            cycle = 0;
            auto = true;
        }
        if is_key_pressed(KeyCode::M) {
            form = match form {
                LeafForm::SmoothPinnate => LeafForm::ToothedPinnate,
                LeafForm::ToothedPinnate => LeafForm::ToothedActinodromous,
                LeafForm::ToothedActinodromous => LeafForm::SmoothPinnate,
            };
            set_leaf_form(form);
            outline_ref = reference_leaf();
            graph = seed_graph();
            sources = throw_darts(
                &graph,
                BirthDistanceSource {
                    ada: BIRTH_DIST_SOURCE,
                },
                BirthDistanceVein {
                    ada: BIRTH_DIST_VEIN,
                },
                DartAttempts { ada: DART_ATTEMPTS },
            );
            t = T0;
            cycle = 0;
            auto = true;
        }
        if is_key_pressed(KeyCode::R) {
            graph = seed_graph();
            sources = throw_darts(
                &graph,
                BirthDistanceSource {
                    ada: BIRTH_DIST_SOURCE,
                },
                BirthDistanceVein {
                    ada: BIRTH_DIST_VEIN,
                },
                DartAttempts { ada: DART_ATTEMPTS },
            );
            t = T0;
            cycle = 0;
            auto = true;
        }
        if auto && t < T_MAX {
            auto_accum += get_frame_time() as f64;
            while auto_accum >= AUTO_STEP_INTERVAL {
                auto_accum -= AUTO_STEP_INTERVAL;
                grow_cycle(&mut graph, &mut sources, &mut t, mode);
                cycle += 1;
            }
        }

        let map_t = model.vertical_map(t, Y0_MAX, N_STEPS);
        let outline = propagate_from_ref(&model, &outline_ref, &map_t0, &map_t, t);

        let fit = fit_world(ref_bbox(T_MAX));

        clear_background(Color::new(0.08, 0.09, 0.11, 1.0));

        draw_polyline(&outline, &fit, leaf_green, 2.5, true);

        for &p in outline.iter().step_by(66) {
            let (sx, sy) = to_screen(&p, &fit);
            draw_circle(sx, sy, 2.5, leaf_green);
        }

        let radii = murray_radii(
            &graph,
            MurrayExponent {
                ada: MURRAY_EXPONENT,
            },
            TipRadius { ada: TIP_RADIUS },
        );
        let max_radius = radii.iter().cloned().fold(0.0f64, f64::max);

        draw_sources(&sources, &fit, source_blue);
        draw_vein_graph(&graph, &radii, &fit, vein_green, node_yellow);

        let mode_name = match mode {
            VenationMode::Open => "open",
            VenationMode::Closed => "closed",
        };
        let form_name = match form {
            LeafForm::SmoothPinnate => "smooth pinnate",
            LeafForm::ToothedPinnate => "toothed pinnate",
            LeafForm::ToothedActinodromous => "toothed actinodromous",
        };
        let loops = (graph.edges.len() + 1).saturating_sub(graph.nodes.len());
        let accel = match mode {
            VenationMode::Closed => "Delaunay",
            VenationMode::Open if graph.nodes.len() >= SPATIAL_INDEX_MIN_NODES => "Delaunay",
            VenationMode::Open => "linear",
        };
        let hud_title = format!(
            "Step 7+8: form = {form_name}   t = {t:.2} / {T_MAX:.1}   cycle = {cycle}   [auto: {}]",
            if auto { "on" } else { "off" }
        );
        draw_text(&hud_title, 20.0, 30.0, 26.0, WHITE);
        let hud_veins = format!(
            "venation = {mode_name} [{accel} NN]  |  nodes = {}  edges = {}  loops = {loops}  n = {MURRAY_EXPONENT:.1}  max r = {max_radius:.1}",
            graph.nodes.len(),
            graph.edges.len()
        );
        draw_text(&hud_veins, 20.0, 56.0, 22.0, node_yellow);
        let phase = if t >= T_MAX {
            "  [mature: blade fully grown]"
        } else {
            ""
        };
        let hud_sources = format!(
            "auxin sources active = {}   (seeded into newly-grown blade each cycle){phase}",
            sources.len()
        );
        draw_text(&hud_sources, 20.0, 82.0, 22.0, source_blue);
        draw_text(
            "[m] leaf form   [o] open/closed   [space] play/pause   [c] one cycle   [r] reset",
            20.0,
            screen_height() - 16.0,
            20.0,
            GRAY,
        );

        next_frame().await
    }
}
