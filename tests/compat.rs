use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-permanova"))
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn parse(table: &str) -> HashMap<String, String> {
    table
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

fn run_ours(dm: &str, grouping: &str, perms: usize, seed: u64) -> String {
    let out = Command::new(ours())
        .arg(golden(dm))
        .args(["-g", &golden(grouping)])
        .args(["-p", &perms.to_string()])
        .args(["--seed", &seed.to_string()])
        .args(["--precision", "9"])
        .output()
        .expect("run rsomics-permanova");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

/// The pseudo-F statistic, sample size, and group count are deterministic and
/// must match scikit-bio to 1e-9. The p-value is a permutation estimate from a
/// different RNG, so it is compared only within Monte-Carlo tolerance.
fn assert_compat(ours: &str, theirs: &str, perms: usize) {
    let a = parse(ours);
    let b = parse(theirs);
    for k in [
        "method name",
        "test statistic name",
        "sample size",
        "number of groups",
        "number of permutations",
    ] {
        assert_eq!(a[k], b[k], "field '{k}' differs");
    }
    let fa: f64 = a["test statistic"].parse().unwrap();
    let fb: f64 = b["test statistic"].parse().unwrap();
    assert!(
        (fa - fb).abs() <= 1e-9 * fb.abs().max(1.0),
        "pseudo-F differs: ours={fa} skbio={fb}"
    );
    if perms > 0 {
        let pa: f64 = a["p-value"].parse().unwrap();
        let pb: f64 = b["p-value"].parse().unwrap();
        // both estimate the same true p with sd ~ sqrt(p(1-p)/perms); allow 5 sd + a floor.
        let sd = (pb * (1.0 - pb) / perms as f64).sqrt();
        let tol = (5.0 * sd).max(3.0 / perms as f64);
        assert!(
            (pa - pb).abs() <= tol,
            "p-value estimate out of Monte-Carlo tolerance: ours={pa} skbio={pb} tol={tol}"
        );
    }
}

fn skbio_python() -> Option<String> {
    for py in ["python3", "python"] {
        let probe = Command::new(py)
            .args(["-c", "import skbio.stats.distance"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py.to_string());
        }
    }
    eprintln!("SKIP live diff: scikit-bio not importable — golden comparison still runs");
    None
}

fn run_oracle(py: &str, dm: &str, grouping: &str, perms: usize, seed: u64) -> String {
    let script = format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"));
    let out = Command::new(py)
        .arg(&script)
        .arg(golden(dm))
        .arg(golden(grouping))
        .arg(perms.to_string())
        .arg(seed.to_string())
        .arg("9")
        .output()
        .expect("run skbio oracle");
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn golden_text(name: &str) -> String {
    std::fs::read_to_string(golden(name)).expect("read golden")
}

const PERMS: usize = 999;
const SEED: u64 = 42;

#[test]
fn statistic_matches_committed_golden() {
    let ours = run_ours("dm_small.tsv", "groups_small.tsv", PERMS, SEED);
    let golden = golden_text("result_small.skbio.tsv");
    assert_compat(&ours, &golden, PERMS);
}

#[test]
fn statistic_matches_committed_golden_three_groups() {
    let ours = run_ours("dm_three.tsv", "groups_three.tsv", PERMS, SEED);
    let golden = golden_text("result_three.skbio.tsv");
    assert_compat(&ours, &golden, PERMS);
}

#[test]
fn statistic_only_matches_golden() {
    let ours = run_ours("dm_small.tsv", "groups_small.tsv", 0, SEED);
    let golden = golden_text("result_small_p0.skbio.tsv");
    assert_compat(&ours, &golden, 0);
}

#[test]
fn live_skbio_diff_small() {
    let Some(py) = skbio_python() else { return };
    let ours = run_ours("dm_small.tsv", "groups_small.tsv", PERMS, SEED);
    let theirs = run_oracle(&py, "dm_small.tsv", "groups_small.tsv", PERMS, SEED);
    assert_compat(&ours, &theirs, PERMS);
}

fn run_ours_expect_failure(dm: &str, grouping: &str) -> (bool, String) {
    let out = Command::new(ours())
        .arg(golden(dm))
        .args(["-g", &golden(grouping)])
        .args(["-p", "0"])
        .output()
        .expect("run rsomics-permanova");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// A distance matrix scikit-bio's `DistanceMatrix` would reject must not yield a
/// confident wrong statistic: the binary bails with a non-zero exit.
#[test]
fn asymmetric_matrix_fails_loud() {
    let (ok, stderr) = run_ours_expect_failure("dm_asymmetric.tsv", "groups_tiny.tsv");
    assert!(
        !ok,
        "asymmetric matrix must exit non-zero; stderr: {stderr}"
    );
    assert!(stderr.contains("symmetric"), "stderr: {stderr}");
}

#[test]
fn nonhollow_matrix_fails_loud() {
    let (ok, stderr) = run_ours_expect_failure("dm_nonhollow.tsv", "groups_tiny.tsv");
    assert!(
        !ok,
        "non-hollow matrix must exit non-zero; stderr: {stderr}"
    );
    assert!(stderr.contains("hollow"), "stderr: {stderr}");
}

#[test]
fn live_skbio_diff_three_groups() {
    let Some(py) = skbio_python() else { return };
    let ours = run_ours("dm_three.tsv", "groups_three.tsv", PERMS, SEED);
    let theirs = run_oracle(&py, "dm_three.tsv", "groups_three.tsv", PERMS, SEED);
    assert_compat(&ours, &theirs, PERMS);
}
