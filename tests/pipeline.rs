//! End-to-end check on synthetic data with two cell types.
//!
//! Real cells are type A or type B, plus a few "true" heterotypic doublets
//! (A + B). With no host-side PCA here, the low-dimensional expression is
//! used directly as the embedding. True doublets should sit among the
//! artificial doublets and get the highest pANN.

use doublet_rs::knn::HnswParams;
use doublet_rs::{
    KnnMethod, call_doublets, k_from_pk, n_doublets_for_pn, pann_from_embedding, simulate_doublets,
};
use ndarray::Array2;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use sprs::{CsMat, TriMat};

const GENES: usize = 20;
const N_A: usize = 300;
const N_B: usize = 300;
const N_TRUE_DOUBLETS: usize = 30;

/// Type A expresses genes 0..10, type B genes 10..20, with noise.
fn profile(rng: &mut ChaCha8Rng, a: f64, b: f64) -> Vec<f64> {
    (0..GENES)
        .map(|g| {
            let base = if g < GENES / 2 { a } else { b };
            (base * 10.0 + rng.random_range(-1.0..1.0)).max(0.0)
        })
        .collect()
}

fn synthetic_cells() -> (CsMat<f64>, Vec<bool>) {
    let mut rng = ChaCha8Rng::seed_from_u64(0);
    let mut cells = Vec::new();
    let mut is_doublet = Vec::new();
    for _ in 0..N_A {
        cells.push(profile(&mut rng, 1.0, 0.0));
        is_doublet.push(false);
    }
    for _ in 0..N_B {
        cells.push(profile(&mut rng, 0.0, 1.0));
        is_doublet.push(false);
    }
    for _ in 0..N_TRUE_DOUBLETS {
        cells.push(profile(&mut rng, 0.5, 0.5));
        is_doublet.push(true);
    }

    let mut tri = TriMat::new((GENES, cells.len()));
    for (c, cell) in cells.iter().enumerate() {
        for (g, &v) in cell.iter().enumerate() {
            if v > 0.0 {
                tri.add_triplet(g, c, v);
            }
        }
    }
    (tri.to_csc(), is_doublet)
}

/// Stacks real cells and artificial doublets into a points × genes array.
fn merged_embedding(real: &CsMat<f64>, doublets: &CsMat<f64>) -> Array2<f64> {
    let n = real.cols() + doublets.cols();
    let mut emb = Array2::zeros((n, GENES));
    for (offset, m) in [(0, real), (real.cols(), doublets)] {
        for (c, col) in m.outer_iterator().enumerate() {
            for (g, &v) in col.iter() {
                emb[[offset + c, g]] = v;
            }
        }
    }
    emb
}

fn run(method: KnnMethod) -> (Vec<f64>, Vec<bool>) {
    let (expr, is_doublet) = synthetic_cells();
    let n_real = expr.cols();

    let n_sim = n_doublets_for_pn(n_real, 0.25).unwrap();
    let doublets = simulate_doublets(expr.view(), n_sim, 1).unwrap();
    let emb = merged_embedding(&expr, &doublets.matrix);

    let k = k_from_pk(emb.nrows(), 0.02);
    let pann = pann_from_embedding(emb.view(), n_real, k, method).unwrap();
    (pann, is_doublet)
}

fn check(pann: &[f64], is_doublet: &[bool]) {
    let mean = |want: bool| {
        let v: Vec<f64> = pann
            .iter()
            .zip(is_doublet)
            .filter(|(_, d)| **d == want)
            .map(|(p, _)| *p)
            .collect();
        v.iter().sum::<f64>() / v.len() as f64
    };
    // Singlets are not near 0: about half the artificial doublets are
    // homotypic (A + A, B + B) and land inside the real clusters.
    let (doublet_mean, singlet_mean) = (mean(true), mean(false));
    assert!(
        doublet_mean > 0.9 && singlet_mean < 0.5,
        "doublet mean pANN {doublet_mean}, singlet mean pANN {singlet_mean}"
    );

    let calls = call_doublets(pann, N_TRUE_DOUBLETS);
    let caught = calls
        .iter()
        .zip(is_doublet)
        .filter(|(c, d)| **c && **d)
        .count();
    assert!(
        caught >= N_TRUE_DOUBLETS * 9 / 10,
        "caught {caught} of {N_TRUE_DOUBLETS}"
    );
}

#[test]
fn exact_pipeline_finds_true_doublets() {
    let (pann, is_doublet) = run(KnnMethod::Exact);
    check(&pann, &is_doublet);
}

#[test]
fn hnsw_pipeline_finds_true_doublets() {
    let (pann, is_doublet) = run(KnnMethod::Hnsw(HnswParams::default()));
    check(&pann, &is_doublet);
}
