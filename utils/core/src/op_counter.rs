// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

//! Operation counters for the STARK-signature verifier (enabled by the `count-ops` feature).
//!
//! All public items in this module are available only when the `count-ops` feature is active.
//!
//! # Semantics
//! * **F-mul** counts every `BaseElement` multiplication that occurs *outside* an extension-field
//!   multiplication.  Multiplications that are performed internally by `QuadExtension::mul` or
//!   `CubeExtension::mul` are suppressed via a thread-local depth guard, so the two counters are
//!   disjoint.
//! * **G-mul** counts every call to `QuadExtension::mul` / `CubeExtension::mul`.
//! * **RPO** counts every call to `Rp64_256::apply_permutation`.
//!
//! # Usage
//! 1. Call [`reset()`] before starting the verifier.
//! 2. Call [`set_block()`] at each phase boundary inside the verifier.
//! 3. After verification completes, call [`report()`] to retrieve the counts.

use core::sync::atomic::{AtomicU64, AtomicU8, Ordering::Relaxed};

// Block indices
pub const BLOCK_TRANSCRIPT: u8 = 0;
pub const BLOCK_DEEP: u8 = 1;
pub const BLOCK_FRI: u8 = 2;
const NUM_BLOCKS: usize = 3;

// Per-block counters  (index 0 = Transcript, 1 = DEEP, 2 = FRI)
static F_MUL: [AtomicU64; NUM_BLOCKS] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
static G_MUL: [AtomicU64; NUM_BLOCKS] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
static RPO:   [AtomicU64; NUM_BLOCKS] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];

/// Index of the currently active block. Values: 0–2 (see `BLOCK_*` constants).
/// Defaults to `BLOCK_TRANSCRIPT` (0).
static CURRENT_BLOCK: AtomicU8 = AtomicU8::new(BLOCK_TRANSCRIPT);

// Thread-local nesting depth: > 0 while executing inside an extension-field multiplication.
// Using a depth counter (rather than a boolean) handles the case where inv() internally calls
// ExtensibleField::mul, which in turn might call BaseElement::mul.
std::thread_local! {
    static G_MUL_DEPTH: std::cell::Cell<u8> = std::cell::Cell::new(0);
}

/// Switch the active accounting block.
#[inline(always)]
pub fn set_block(block: u8) {
    CURRENT_BLOCK.store(block, Relaxed);
}

/// Mark entry into an extension-field multiplication; suppresses F-mul counting until the
/// matching [`exit_g_mul`] call.
#[inline(always)]
pub fn enter_g_mul() {
    G_MUL_DEPTH.with(|d| d.set(d.get().saturating_add(1)));
}

/// Mark exit from an extension-field multiplication.
#[inline(always)]
pub fn exit_g_mul() {
    G_MUL_DEPTH.with(|d| d.set(d.get().saturating_sub(1)));
}

/// Increment the Goldilocks (base-field) multiplication counter by 1.
/// No-ops when called from within an extension-field multiplication.
#[inline(always)]
pub fn inc_f_mul() {
    if G_MUL_DEPTH.with(|d| d.get()) > 0 {
        return;
    }
    let b = CURRENT_BLOCK.load(Relaxed) as usize;
    F_MUL[b].fetch_add(1, Relaxed);
}

/// Increment the extension-field multiplication counter by 1.
#[inline(always)]
pub fn inc_g_mul() {
    let b = CURRENT_BLOCK.load(Relaxed) as usize;
    G_MUL[b].fetch_add(1, Relaxed);
}

/// Increment the RPO permutation call counter by 1.
#[inline(always)]
pub fn inc_rpo() {
    let b = CURRENT_BLOCK.load(Relaxed) as usize;
    RPO[b].fetch_add(1, Relaxed);
}

/// Zero all counters and reset the active block to `BLOCK_TRANSCRIPT`.
pub fn reset() {
    for i in 0..NUM_BLOCKS {
        F_MUL[i].store(0, Relaxed);
        G_MUL[i].store(0, Relaxed);
        RPO[i].store(0, Relaxed);
    }
    CURRENT_BLOCK.store(BLOCK_TRANSCRIPT, Relaxed);
}

// REPORT TYPES
// ================================================================================================

/// Counts for a single verifier block.
#[derive(Debug, Clone, Copy, Default)]
pub struct BlockCounts {
    /// Base-field multiplications *not* occurring inside an extension-field multiplication.
    pub f_mul: u64,
    /// Extension-field multiplications.
    pub g_mul: u64,
    /// RPO permutation calls.
    pub rpo:   u64,
}

impl BlockCounts {
    fn read(idx: usize) -> Self {
        Self {
            f_mul: F_MUL[idx].load(Relaxed),
            g_mul: G_MUL[idx].load(Relaxed),
            rpo:   RPO[idx].load(Relaxed),
        }
    }

    fn total(self, other: Self) -> Self {
        Self {
            f_mul: self.f_mul + other.f_mul,
            g_mul: self.g_mul + other.g_mul,
            rpo:   self.rpo   + other.rpo,
        }
    }

    /// F-mul cost induced by the G-mul operations, given the Karatsuba coefficient `c_e`
    /// (the number of base-field multiplications per extension-field multiplication).
    pub fn g_mul_induced_f_mul(&self, c_e: u64) -> u64 {
        self.g_mul * c_e
    }
}

/// Per-block and total operation counts for one verifier execution.
#[derive(Debug, Clone, Copy)]
pub struct CountReport {
    pub transcript: BlockCounts,
    pub deep:       BlockCounts,
    pub fri:        BlockCounts,
    pub total:      BlockCounts,
}

/// Snapshot the current counter values and return a [`CountReport`].
pub fn report() -> CountReport {
    let transcript = BlockCounts::read(BLOCK_TRANSCRIPT as usize);
    let deep       = BlockCounts::read(BLOCK_DEEP as usize);
    let fri        = BlockCounts::read(BLOCK_FRI as usize);
    let total      = transcript.total(deep).total(fri);
    CountReport { transcript, deep, fri, total }
}
