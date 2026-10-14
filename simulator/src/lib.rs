//! CPU simulator for the VIRAL mean-field QAOA FPGA design.
//!
//! - `hw`, `cordic`: bit-exact model of VIRAL_n_spin (rank-1 coupling)
//! - `model`: floating-point model with rank-1 or dense coupling
//! - `problem`, `baseline`: number partitioning and SK spin glasses, and
//!   classical methods to compare against
//! - `run`, `server`: runs for the web page in ../web
//! - `bench`: CPU timings against the FPGA estimates

pub mod baseline;
pub mod bench;
pub mod cordic;
pub mod fft;
pub mod hexfile;
pub mod hw;
pub mod json;
pub mod model;
pub mod par;
pub mod problem;
pub mod rng;
pub mod run;
pub mod server;
pub mod study;
