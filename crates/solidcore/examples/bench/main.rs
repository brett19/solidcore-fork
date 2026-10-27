//! `cargo run --profile difftest -p solidcore --example bench -- [filter] [--reps N] [--json <path>]`:
//! end-to-end performance cases built in pure Rust through the public API, for
//! profiling and for tracking speed over time.
//!
//! Each case is a recipe (no input files): a common modelling pattern (hole plates,
//! patterned cuts, filleted and chamfered boxes, an enclosure) or a known slow one
//! (coplanar cuts of a Menger sponge, curved unions of spheres and tori, fillets on a
//! boolean result). Every case is built twice in fresh kernels: the first run (or the
//! fastest of `--reps` runs) is timed per stage, and the second must give the same face
//! count and the same volume and area bits (`det`), so a scheduling- or address-dependent
//! result shows up here. A case that fails reports its error instead of a time; a failure
//! is a result too.
//!
//! Stages: `build` (primitives and transforms), `boolean`, `fillet` (fillets and
//! chamfers), `volume` and `area` (`brep_volume`, `brep_area` at accuracy 1e-6).
//!
//! `--json` writes one object per case:
//! `{ "name", "ok", "error", "ms", "stages": { "<stage>": ms }, "faces", "solid", "det", "volume", "area" }`.

use solidcore::blend::fillet_executive::api::{chamfer_fillet, fillet_edges};
use solidcore::boolean::api_brep::{boolean_difference, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::mass::integrator::properties::{brep_area, brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::{
    create_box, create_cylinder, create_sphere, create_torus,
};
use solidcore::topology::{edge, vertex, EdgeId};
use solidcore::{Brep, Kernel, Result};
use std::f64::consts::PI;
use std::fmt::Write as _;
use std::time::{Duration, Instant};

/// Time spent per stage, in first-use order.
#[derive(Default)]
struct Stages(Vec<(&'static str, Duration)>);

impl Stages {
    /// Runs `f`, adding its time to `stage`.
    fn time<T>(&mut self, stage: &'static str, f: impl FnOnce() -> Result<T>) -> Result<T> {
        let t = Instant::now();
        let r = f();
        let d = t.elapsed();
        match self.0.iter_mut().find(|s| s.0 == stage) {
            Some(s) => s.1 += d,
            None => self.0.push((stage, d)),
        }
        r
    }
}

/// A recipe: builds its body in `k`.
type Recipe = fn(&mut Kernel, &mut Stages) -> Result<Brep>;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Rotates `b` by `angle` about the z axis.
fn rotate_z(k: &mut Kernel, b: &mut Brep, angle: f64) -> Result<()> {
    let placement = Axis2Placement::new(
        p(0.0, 0.0, 0.0),
        Vec3::new(angle.cos(), angle.sin(), 0.0),
        Vec3::new(-angle.sin(), angle.cos(), 0.0),
    );
    solidcore::topology::brep::p5::transform(k, b, &placement, None, false)
}

/// The edges of `b` whose end points satisfy `keep`.
fn edges_where(b: &Brep, keep: impl Fn(Point3, Point3) -> bool) -> Result<Vec<EdgeId>> {
    let mut out = Vec::new();
    for e in b.edges().ids() {
        if let (Some(s), Some(t)) = edge::start_end_vertices(b, e)? {
            if keep(vertex::point(b, s)?, vertex::point(b, t)?) {
                out.push(e);
            }
        }
    }
    Ok(out)
}

/// Whether the edge from `a` to `c` is parallel to z.
fn vertical(a: Point3, c: Point3) -> bool {
    (a.x - c.x).abs() < 1e-9 && (a.y - c.y).abs() < 1e-9
}

/// A 100 x 100 x 10 plate with `n` x `n` r 2 through holes, cut one at a time.
fn plate_holes(k: &mut Kernel, st: &mut Stages, n: usize) -> Result<Brep> {
    let mut plate = st.time("build", || {
        create_box(k, p(0.0, 0.0, 0.0), 100.0, 100.0, 10.0)
    })?;
    let pitch = 100.0 / n as f64;
    for i in 0..n {
        for j in 0..n {
            let centre = (pitch * (i as f64 + 0.5), pitch * (j as f64 + 0.5));
            let tool = st.time("build", || {
                create_cylinder(k, p(centre.0, centre.1, -1.0), 2.0, 12.0)
            })?;
            plate = st.time("boolean", || boolean_difference(k, plate, tool))?;
        }
    }
    Ok(plate)
}

/// An r 50 disc minus `z` rectangular tooth gaps, rotated about the axis (a patterned cut).
fn gear(k: &mut Kernel, st: &mut Stages, z: usize) -> Result<Brep> {
    let mut disc = st.time("build", || create_cylinder(k, p(0.0, 0.0, 0.0), 50.0, 10.0))?;
    let gap = 2.0 * PI * 45.0 / z as f64 * 0.45;
    for i in 0..z {
        let angle = 2.0 * PI * i as f64 / z as f64 + 0.1234;
        let mut tool = st.time("build", || {
            create_box(k, p(45.0, -gap / 2.0, -1.0), 10.0, gap, 12.0)
        })?;
        st.time("build", || rotate_z(k, &mut tool, angle))?;
        disc = st.time("boolean", || boolean_difference(k, disc, tool))?;
    }
    Ok(disc)
}

/// A 40 x 30 x 20 box with its 4 vertical edges (or all 12) filleted r 3.
fn filleted_box(k: &mut Kernel, st: &mut Stages, all: bool) -> Result<Brep> {
    let mut b = st.time("build", || {
        create_box(k, p(0.0, 0.0, 0.0), 40.0, 30.0, 20.0)
    })?;
    let edges = if all {
        b.edges().ids().collect()
    } else {
        edges_where(&b, vertical)?
    };
    st.time("fillet", || fillet_edges(k, &mut b, &edges, 3.0, 1, 1, 1.0))?;
    Ok(b)
}

/// The same box with every edge chamfered 2 (`SmApiChamferFillet`).
fn chamfer_box(k: &mut Kernel, st: &mut Stages) -> Result<Brep> {
    let mut b = st.time("build", || {
        create_box(k, p(0.0, 0.0, 0.0), 40.0, 30.0, 20.0)
    })?;
    st.time("fillet", || chamfer_fillet(k, &mut b, 2.0))?;
    Ok(b)
}

/// An enclosure body: a box with its vertical edges filleted r 6, four corner bosses
/// joined and screw holes cut into them, then `vents` x `vents` vent holes.
fn enclosure(k: &mut Kernel, st: &mut Stages, vents: usize) -> Result<Brep> {
    let mut b = st.time("build", || {
        create_box(k, p(0.0, 0.0, 0.0), 80.0, 60.0, 30.0)
    })?;
    let edges = edges_where(&b, vertical)?;
    st.time("fillet", || fillet_edges(k, &mut b, &edges, 6.0, 1, 1, 1.0))?;
    for (x, y) in [(12.0, 12.0), (68.0, 12.0), (12.0, 48.0), (68.0, 48.0)] {
        let boss = st.time("build", || create_cylinder(k, p(x, y, 30.0), 4.0, 5.0))?;
        b = st.time("boolean", || boolean_union(k, b, boss))?;
        let hole = st.time("build", || create_cylinder(k, p(x, y, 20.0), 1.5, 16.0))?;
        b = st.time("boolean", || boolean_difference(k, b, hole))?;
    }
    for i in 0..vents {
        for j in 0..vents {
            let x = 25.0 + 30.0 * (i as f64 + 0.5) / vents as f64;
            let y = 15.0 + 30.0 * (j as f64 + 0.5) / vents as f64;
            let vent = st.time("build", || create_cylinder(k, p(x, y, -1.0), 0.8, 3.0))?;
            b = st.time("boolean", || boolean_difference(k, b, vent))?;
        }
    }
    Ok(b)
}

/// A Menger sponge of level 1 or 2 cut from a 27 cube by through-bars along each axis:
/// every cut meets earlier ones in coplanar faces.
fn menger(k: &mut Kernel, st: &mut Stages, level: u32) -> Result<Brep> {
    let s = 27.0;
    let mut b = st.time("build", || create_box(k, p(0.0, 0.0, 0.0), s, s, s))?;
    let mut bars = vec![(9.0, 9.0, 9.0)];
    if level >= 2 {
        for i in 0..3 {
            for j in 0..3 {
                if (i, j) != (1, 1) {
                    bars.push((f64::from(i) * 9.0 + 3.0, f64::from(j) * 9.0 + 3.0, 3.0));
                }
            }
        }
    }
    for axis in 0..3 {
        for &(a, c, w) in &bars {
            let (origin, size) = match axis {
                0 => (p(-1.0, a, c), (s + 2.0, w, w)),
                1 => (p(a, -1.0, c), (w, s + 2.0, w)),
                _ => (p(a, c, -1.0), (w, w, s + 2.0)),
            };
            let tool = st.time("build", || create_box(k, origin, size.0, size.1, size.2))?;
            b = st.time("boolean", || boolean_difference(k, b, tool))?;
        }
    }
    Ok(b)
}

/// `n` x `n` overlapping r 1 spheres at spacing 1.5, joined one at a time.
fn sphere_grid(k: &mut Kernel, st: &mut Stages, n: usize) -> Result<Brep> {
    let mut b = st.time("build", || create_sphere(k, p(0.0, 0.0, 0.0), 1.0))?;
    for i in 0..n {
        for j in 0..n {
            if (i, j) != (0, 0) {
                let c = p(1.5 * i as f64, 1.5 * j as f64, 0.0);
                let s = st.time("build", || create_sphere(k, c, 1.0))?;
                b = st.time("boolean", || boolean_union(k, b, s))?;
            }
        }
    }
    Ok(b)
}

/// `n` tori (R 5, r 1.5, centred 4 off the z axis) rotated about z and joined: torus/torus
/// intersections.
fn torus_ring(k: &mut Kernel, st: &mut Stages, n: usize) -> Result<Brep> {
    let mut make = |k: &mut Kernel, i: usize| -> Result<Brep> {
        let mut t = st.time("build", || create_torus(k, p(4.0, 0.0, 0.0), 5.0, 1.5))?;
        let angle = 2.0 * PI * i as f64 / n as f64 + 0.05;
        st.time("build", || rotate_z(k, &mut t, angle))?;
        Ok(t)
    };
    let mut b = make(k, 0)?;
    let mut rest = Vec::new();
    for i in 1..n {
        rest.push(make(k, i)?);
    }
    for t in rest {
        b = st.time("boolean", || boolean_union(k, b, t))?;
    }
    Ok(b)
}

/// The 4 x 4 hole plate with the top edge of every hole filleted r 0.5.
fn fillet_hole_tops(k: &mut Kernel, st: &mut Stages) -> Result<Brep> {
    let mut b = plate_holes(k, st, 4)?;
    let inside = |q: Point3| {
        (q.z - 10.0).abs() < 1e-9 && q.x > 0.0 && q.x < 100.0 && q.y > 0.0 && q.y < 100.0
    };
    let edges = edges_where(&b, |a, c| inside(a) && inside(c))?;
    st.time("fillet", || fillet_edges(k, &mut b, &edges, 0.5, 1, 1, 1.0))?;
    Ok(b)
}

/// Every case, by name.
fn cases() -> Vec<(&'static str, Recipe)> {
    vec![
        ("plate_holes_4", |k, st| plate_holes(k, st, 4)),
        ("plate_holes_8", |k, st| plate_holes(k, st, 8)),
        ("plate_holes_12", |k, st| plate_holes(k, st, 12)),
        ("gear_20", |k, st| gear(k, st, 20)),
        ("gear_60", |k, st| gear(k, st, 60)),
        ("fillet_box_4", |k, st| filleted_box(k, st, false)),
        ("fillet_box_12", |k, st| filleted_box(k, st, true)),
        ("chamfer_box", chamfer_box),
        ("enclosure_3", |k, st| enclosure(k, st, 3)),
        ("enclosure_6", |k, st| enclosure(k, st, 6)),
        ("fillet_holes_4", fillet_hole_tops),
        ("menger_1", |k, st| menger(k, st, 1)),
        ("menger_2", |k, st| menger(k, st, 2)),
        ("sphere_grid_3", |k, st| sphere_grid(k, st, 3)),
        ("sphere_grid_5", |k, st| sphere_grid(k, st, 5)),
        ("torus_ring_4", |k, st| torus_ring(k, st, 4)),
        ("torus_ring_8", |k, st| torus_ring(k, st, 8)),
    ]
}

/// What a run of a case measured.
struct Measured {
    stages: Stages,
    faces: usize,
    solid: bool,
    volume: Option<f64>,
    area: Option<f64>,
}

/// Builds `recipe` in a fresh kernel, then measures the body.
fn run(recipe: Recipe) -> Result<Measured> {
    let mut k = Kernel::default();
    let mut stages = Stages::default();
    let mut brep = recipe(&mut k, &mut stages)?;
    let faces = solidcore::topology::brep::p1::faces(&brep, None)?.len();
    let solid = solidcore::topology::brep::p2::is_manifold_solid(&brep)?;
    let store = AttributeStore::new();
    let key_of = |_| ObjectKey(0);
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &key_of,
    };
    let b = &mut brep;
    let volume = stages
        .time("volume", || brep_volume(&k, b, 1e-6, &attrs))
        .ok();
    let area = stages.time("area", || brep_area(&k, b, 1e-6)).ok();
    Ok(Measured {
        stages,
        faces,
        solid,
        volume,
        area,
    })
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1e3
}

fn bits(x: Option<f64>) -> Option<u64> {
    x.map(f64::to_bits)
}

fn json_number(x: Option<f64>) -> String {
    x.filter(|v| v.is_finite())
        .map_or_else(|| "null".into(), |v| format!("{v:?}"))
}

fn main() -> std::result::Result<(), String> {
    let mut filter = String::new();
    let mut reps = 1_usize;
    let mut json_path = None;
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--reps" => {
                reps = args
                    .next()
                    .and_then(|n| n.parse().ok())
                    .filter(|&n| n > 0)
                    .ok_or("--reps needs a positive count")?;
            }
            "--json" => json_path = Some(args.next().ok_or("--json needs a path")?),
            _ => filter = a,
        }
    }
    let mut json = Vec::new();
    for (name, recipe) in cases() {
        if !name.contains(&filter) {
            continue;
        }
        // the fastest of `reps` runs, then one more for the determinism check
        let mut best: Option<(Duration, Result<Measured>)> = None;
        for _ in 0..reps {
            let t = Instant::now();
            let r = run(recipe);
            let d = t.elapsed();
            if best.as_ref().is_none_or(|(b, _)| d < *b) {
                best = Some((d, r));
            }
        }
        let Some((wall, first)) = best else { continue };
        let second = run(recipe);
        let mut entry = format!("{{\"name\": \"{name}\", \"ms\": {:.1}", ms(wall));
        match (first, second) {
            (Ok(a), Ok(b)) => {
                let det = a.faces == b.faces
                    && bits(a.volume) == bits(b.volume)
                    && bits(a.area) == bits(b.area);
                let stages: Vec<String> = a
                    .stages
                    .0
                    .iter()
                    .map(|(s, d)| format!("{s} {:.0}", ms(*d)))
                    .collect();
                println!(
                    "{name:<16} {:>8.0} ms  faces {:<4} solid {:<5} det {:<5} volume {} area {}  [{}]",
                    ms(wall),
                    a.faces,
                    a.solid,
                    det,
                    json_number(a.volume),
                    json_number(a.area),
                    stages.join(", "),
                );
                let stage_json: Vec<String> = a
                    .stages
                    .0
                    .iter()
                    .map(|(s, d)| format!("\"{s}\": {:.1}", ms(*d)))
                    .collect();
                let _ = write!(
                    entry,
                    ", \"ok\": true, \"error\": null, \"stages\": {{{}}}, \"faces\": {}, \"solid\": {}, \"det\": {det}, \"volume\": {}, \"area\": {}}}",
                    stage_json.join(", "),
                    a.faces,
                    a.solid,
                    json_number(a.volume),
                    json_number(a.area),
                );
            }
            (Err(e), _) | (_, Err(e)) => {
                println!("{name:<16} {:>8.0} ms  FAILED: {e}", ms(wall));
                let msg = e.to_string().replace('\\', "\\\\").replace('"', "\\\"");
                let _ = write!(
                    entry,
                    ", \"ok\": false, \"error\": \"{msg}\", \"stages\": {{}}, \"faces\": null, \"solid\": null, \"det\": null, \"volume\": null, \"area\": null}}"
                );
            }
        }
        json.push(entry);
    }
    if let Some(path) = json_path {
        std::fs::write(&path, format!("[\n{}\n]\n", json.join(",\n")))
            .map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}
