extern crate bigsi_rs;
extern crate bincode;
extern crate bv;
extern crate xxh3;
extern crate flate2;
extern crate fnv;
extern crate probability;
extern crate rayon;
extern crate serde;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate itertools;

pub mod build;

pub mod read_id_mt_pe;

pub mod read_id_batch;

pub mod perfect_search;

pub mod batch_search_pe;

pub mod kmer;

pub mod reports;

pub mod read_filter;

pub mod seq;
