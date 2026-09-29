//! End-to-end check on synthetic data with two cell types.
//!
//! Real cells are type A or type B, plus a few "true" heterotypic doublets
//! (A + B). With no host-side PCA here, the low-dimensional expression is
//! used directly as the embedding. True doublets should sit among the
//! artificial doublets and get the highest pANN.

use doublet_rs::{
    call_doublets, k_from_pk, n_doublets_for_pn, pann_from_embedding, simulate_doublets, CscMatrix,
    Embedding, KnnMethod,
};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

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

/// Dense cells (each a GENES-long vector) and which are true doublets.
fn synthetic_cells() -> (Vec<Vec<f64>>, Vec<bool>) {
    let mut rng = ChaCha8Rng::seed_from_u64(0);
    let mut cells = Vec::new();
    let mut is_doublet = Vec::new();
    for (n, a, b, doublet) in [
        (N_A, 1.0, 0.0, false),
        (N_B, 0.0, 1.0, false),
        (N_TRUE_DOUBLETS, 0.5, 0.5, true),
    ] {
        for _ in 0..n {
            cells.push(profile(&mut rng, a, b));
            is_doublet.push(doublet);
        }
    }
    (cells, is_doublet)
}

fn to_csc(cells: &[Vec<f64>]) -> CscMatrix<usize> {
    let triplets = cells
        .iter()
        .enumerate()
        .flat_map(|(c, cell)| {
            cell.iter()
                .enumerate()
                .filter(|(_, v)| **v > 0.0)
                .map(move |(g, &v)| (g, c, v))
        })
        .collect();
    CscMatrix::from_triplets(GENES, cells.len(), triplets).unwrap()
}

fn run(method: KnnMethod) -> (Vec<f64>, Vec<bool>) {
    let (cells, is_doublet) = synthetic_cells();
    let expr = to_csc(&cells);
    let n_real = cells.len();

    let n_sim = n_doublets_for_pn(n_real, 0.25).unwrap();
    let doublets = simulate_doublets(expr.view(), n_sim, 1).unwrap();

    // Row-major merged embedding: real cells, then artificial doublets.
    let mut merged: Vec<f64> = cells.concat();
    let dv = doublets.matrix.view();
    for c in 0..n_sim {
        merged.extend((0..GENES).map(|g| dv.get(g, c)));
    }
    let n_points = n_real + n_sim;
    let emb = Embedding::new(&merged, n_points, GENES).unwrap();

    let k = k_from_pk(n_points, 0.02);
    let pann = pann_from_embedding(emb, n_real, k, method).unwrap();
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

#[cfg(feature = "hnsw")]
#[test]
fn hnsw_pipeline_finds_true_doublets() {
    let (pann, is_doublet) = run(KnnMethod::Hnsw(doublet_rs::knn::HnswParams::default()));
    check(&pann, &is_doublet);
}
