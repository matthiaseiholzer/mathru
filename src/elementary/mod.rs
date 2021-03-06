//! Elementary functions
//!
//! Fore more information:
//! <a href="https://en.wikipedia.org/wiki/Elementary_function">https://en.wikipedia.org/wiki/Elementary_function</a>

mod power;

mod exponential;

mod trigonometry;

mod hyperbolic;

pub use self::{
    exponential::Exponential, hyperbolic::Hyperbolic, power::Power, trigonometry::Trigonometry,
};
