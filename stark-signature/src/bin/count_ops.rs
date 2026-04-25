//! Measurement binary that counts the operations performed by the STARK-signature verifier
//! and prints the data needed to fill Table~\ref{tab:verification-cost} in performance.tex.
//!
//! Run with:
//!   cargo run --bin count_ops --features count-ops
//!
//! The binary measures one execution of SigVer for each of the two parameter regimes
//! (UDR and LDR) with the BCS hash instantiated as RPO (Rp64_256), and reports:
//!   - Per-block counts: Transcript hashing / DEEP algebraic check / FRI query phase
//!   - Totals
//!
//! F-mul counts only base-field multiplications that occur *outside* an extension-field
//! multiplication (the two counters are disjoint).  The induced F-mul cost of the G-muls
//! is displayed separately as c_e × G-mul, where c_e is the number of base-field
//! multiplications per extension-field multiplication under the Karatsuba formulae used:
//!   e=2 (QuadExtension): c_e = 3
//!   e=3 (CubeExtension): c_e = 6

use air::{FieldExtension, ProofOptions};
use crypto::hashers::Rp64_256;
use math::{fields::f64::BaseElement, FieldElement, StarkField};
use rand::Rng;
use stark_signature::{hash, RpoSignature};
use utils::op_counter::{self, CountReport};

// ================================================================================================
// PARAMETER SETS
// ================================================================================================

/// Unique decoding regime: n_FRI = 126, e = 2 (quadratic extension), t0 = 2.
fn udr_options() -> ProofOptions {
    ProofOptions::new(126, 8, 0, FieldExtension::Quadratic, 2, 255, true)
}

/// List decoding regime: n_FRI = 85, e = 3 (cubic extension), t0 = 8.
fn ldr_options() -> ProofOptions {
    ProofOptions::new(85, 8, 0, FieldExtension::Cubic, 8, 255, true)
}

// ================================================================================================
// MAIN
// ================================================================================================

fn main() {
    let mut rng = rand::thread_rng();

    // Generate a random secret key and message once; reuse for both measurements.
    let sk: [BaseElement; 4] = random_array(&mut rng);
    let msg: [BaseElement; 4] = random_array(&mut rng);
    let pk = hash(sk);

    println!("STARK-signature verifier operation counts");
    println!("BCS hash: RPO (Rp64_256)");
    println!("F-mul: base-field muls outside G-mul; G-mul: extension-field muls; c_e*G-mul: induced F-mul cost");
    println!();

    // c_e = Karatsuba base-field muls per G-mul: 3 for e=2, 6 for e=3.
    let udr = measure("Unique decoding regime (n_FRI=126, e=2, t0=2)", sk, msg, pk, udr_options(), 3);
    let ldr = measure("List decoding regime   (n_FRI=85,  e=3, t0=8)", sk, msg, pk, ldr_options(), 6);

    println!();
    println!("=== LaTeX table rows (copy into tab:verification-cost) ===");
    println!();
    println!("% Unique decoding regime (n_FRI=126, e=2, t0=2, c_e=3)");
    print_latex_rows(&udr, 3);
    println!();
    println!("% List decoding regime (n_FRI=85, e=3, t0=8, c_e=6)");
    print_latex_rows(&ldr, 6);
}

// ================================================================================================
// MEASUREMENT
// ================================================================================================

fn measure(
    label: &str,
    sk: [BaseElement; 4],
    msg: [BaseElement; 4],
    pk: [BaseElement; 4],
    options: ProofOptions,
    c_e: u64,
) -> CountReport {
    let sig: RpoSignature<Rp64_256> = RpoSignature::new(options.clone());

    // Sign (prover) — not measured.
    let proof = sig.sign(sk, msg);

    // Reset counters, then verify (verifier) — measured.
    op_counter::reset();
    let result = sig.verify(pk, msg, proof);
    let report = op_counter::report();

    assert!(result.is_ok(), "Verification failed — counts may be incomplete");

    println!("--- {} ---", label);
    println!("{:<30} {:>12} {:>10} {:>12} {:>8}",
        "Block", "F-mul", "G-mul", "c_e*G-mul", "RPO");
    println!("{}", "-".repeat(78));
    print_block("Transcript hashing",   &report.transcript, c_e);
    print_block("DEEP algebraic check", &report.deep,       c_e);
    print_block("FRI query phase",      &report.fri,        c_e);
    println!("{}", "-".repeat(78));
    print_block("Total",                &report.total,      c_e);
    println!();

    report
}

fn print_block(name: &str, b: &utils::op_counter::BlockCounts, c_e: u64) {
    println!("{:<30} {:>12} {:>10} {:>12} {:>8}",
        name,
        b.f_mul,
        b.g_mul,
        b.g_mul_induced_f_mul(c_e),
        b.rpo,
    );
}

fn print_latex_rows(report: &CountReport, c_e: u64) {
    let b = &report.transcript;
    println!("Transcript hashing   & {:>12} & {:>8} & {:>10} & {:>5} \\\\",
        b.f_mul, b.g_mul, b.g_mul_induced_f_mul(c_e), b.rpo);
    let b = &report.deep;
    println!("DEEP algebraic check & {:>12} & {:>8} & {:>10} & {:>5} \\\\",
        b.f_mul, b.g_mul, b.g_mul_induced_f_mul(c_e), b.rpo);
    let b = &report.fri;
    println!("FRI query phase      & {:>12} & {:>8} & {:>10} & {:>5} \\\\",
        b.f_mul, b.g_mul, b.g_mul_induced_f_mul(c_e), b.rpo);
    let b = &report.total;
    println!("\\emph{{Total}}       & {:>12} & {:>8} & {:>10} & {:>5} \\\\",
        b.f_mul, b.g_mul, b.g_mul_induced_f_mul(c_e), b.rpo);
}

// ================================================================================================
// HELPERS
// ================================================================================================

fn random_array<R: Rng>(rng: &mut R) -> [BaseElement; 4] {
    let mut out = [BaseElement::ZERO; 4];
    for x in out.iter_mut() {
        let raw: u64 = rng.gen_range(0..BaseElement::MODULUS);
        *x = BaseElement::new(raw);
    }
    out
}
