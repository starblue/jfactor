//! A library for factoring integers and computing various related functions.

#![warn(missing_docs)]

#![cfg_attr(all(feature = "bench", test), feature(test))]

#[cfg(all(feature = "bench", test))]
extern crate test;

extern crate num;


mod factor;
mod divisors;
mod multiplicative_functions;

pub use factor::*;
