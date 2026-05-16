//! M1 acceptance: the canonical schema (docs/new/simulator_schema.json)
//! and the embedded copy in src/synth/simulator.schema.json must stay
//! byte-for-byte identical.

#[test]
fn embedded_schema_matches_canonical() {
    let canonical = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("docs")
            .join("new")
            .join("simulator_schema.json"),
    )
    .expect("read docs/new/simulator_schema.json");
    let embedded = kitehor::synth::CANONICAL_SCHEMA;
    assert_eq!(
        canonical, embedded,
        "schema drift: docs/new/simulator_schema.json differs from src/synth/simulator.schema.json — \
         copy the canonical file into src/synth/ to resync"
    );
}

#[test]
fn embedded_schema_is_valid_json() {
    let v: serde_json::Value = serde_json::from_str(kitehor::synth::CANONICAL_SCHEMA)
        .expect("embedded schema must parse as JSON");
    assert!(v.get("$schema").is_some(), "missing $schema field");
    assert!(v.get("$defs").is_some(), "missing $defs section");
}
