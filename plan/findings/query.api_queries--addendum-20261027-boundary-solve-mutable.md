Updates: plan/units/query.api_queries.toml (notes: "Partial since the v1.0.1 pin move")

# query.api_queries: v1.0.1 deferral resolved (boundary solve on the mutable B-rep)

Unit: query.api_queries. Inputs: the 16 `relationship_*` op7/op8 cases deferred since the v1.0.1 pin.

v1.0.1 `sm_CountContainedSamples` (SmApiQueries.cpp:139) calls `Point3DClassify(..., TRUE)`, turning
`bDoBoundaryIntersections` on. `BrepQuery::classify` solved a shared target (`BrepClassifier`) on the
read-only `BrepObjects`; the surface cache that solve built had no face binding, so the cases deferred at
`SurfaceCache::face` (surface_cache/p1.rs, `Unported topology.face_p2`). That site is correct for an unbound
cache; the defect was the missing binding. Upstream's solve always has the face (`smsurf_GetSurfaceCache`).

Fix: the boundary solve runs on `target.mutable()` with `solve_mutable_brep_point`, as the ray fires already
do. Numerics and order unchanged. Unit test `classification...`/`classifier` assertion `has_copied` now
expects a copy, since the boundary solve copies first. Result: 1203/1203 pass, 0 deferred; unit is green.
