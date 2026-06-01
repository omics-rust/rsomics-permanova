use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

fn bench_permanova(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-permanova");
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let dm = dir.join("dm_small.tsv");
    let grp = dir.join("groups_small.tsv");
    c.bench_function("rsomics-permanova small 999perm", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(&dm)
                .args(["-g", grp.to_str().unwrap()])
                .args(["-p", "999", "--seed", "42", "-t", "1"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_permanova);
criterion_main!(benches);
