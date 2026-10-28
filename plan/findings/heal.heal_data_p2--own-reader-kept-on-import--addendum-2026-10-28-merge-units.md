# boolean.merge_p2, boolean.merge_p5: operands now built by the Rust kernel (2026-10-28)

Updates: plan/findings/heal.heal_data_p2--own-reader-kept-on-import.md

Scope: the two merge units only (`boolean.merge_p2`, `boolean.merge_p5`); the
`heal.heal_data_p2` part of that finding is handled separately.

The earlier finding called the merge operands "recorded graph fixtures" needing a
fixture rewrite. They are not recorded graphs. Every `graph` key in
`fixtures/boolean.merge_p2/cases.json` and `fixtures/boolean.merge_p5/*.json` is already a
construction recipe (`Spec::Graph { brep }`, a `topology.vertexuse` `BrepIn`: an `.smb`
file, a wire edge from a curve, a face from a surface, or the `prog_test` wire network).
The shim (`sc_merge_fixture.h` `operand`) builds each one natively with the shared fixture
`build`. Only the Rust side went through C++: it called the shim's mode `-1` export and
imported that graph together with the vertex use lists.

Now the Rust side builds the same recipe itself (`boolean_merge_p5::own_operand`). An empty
B-rep and a shell vertex are built with the Rust constructors. Every other kind goes
through `crate::own_builder::build_any`, which uses the Rust `.smb` reader or the Rust
constructors. The result is renumbered into the shim's `native_graph` order (the edge list,
then the vertex list, then the `sc_topology_fixture.h` walk) with
`own_builder::renumber_native_lists`. The mate defaults and the generator's edge counts
now come from the Rust build. The generator first applies the shim's pure input guard
(`solidcore_ffi::units::boolean_merge_p5::check_operand`), so it still skips operands the
shim refuses, as the C++ export used to make it do.

The fixtures did not change, and every case id is kept. Before and after on the same base:
`boolean.merge_p2` 45/45 and `boolean.merge_p5` 140/140 pass. The per-case Rust output
digests (`DIFFTEST_RUST_DIGEST`) are bit-identical before and after for all 185 cases, so
no coverage was lost and no mismatch appeared.
