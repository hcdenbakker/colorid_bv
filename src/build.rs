//extern crate bit_vec;
//extern crate rayon;

use super::kmer;
use bigsi_rs;
use fnv;
use rayon::prelude::*;
use std;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::sync::Mutex;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct Parameters {
    pub k_size: usize,
    pub m_size: usize,
}

pub fn tab_to_map(filename: String) -> fnv::FnvHashMap<std::string::String, Vec<String>> {
    let mut map = fnv::FnvHashMap::default();
    let f = File::open(filename).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        if v.len() == 2 {
            map.insert(String::from(v[0]), vec![String::from(v[1])]);
        } else {
            map.insert(
                String::from(v[0]),
                vec![String::from(v[1]), String::from(v[2])],
            );
        }
    }
    map
}

#[allow(unused_assignments)]
pub fn build_multi_mini_bigvec(
    map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: u64,
    k_size: usize,
    m: usize,
    t: usize,
    _quality: u8,
    _cutoff: isize,
    batch: usize,
) -> (
    bigsi_rs::Bigsi,
    fnv::FnvHashMap<usize, String>,
    fnv::FnvHashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    rayon::ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let map_length = map.len();
    let mut bigsi = bigsi_rs::Bigsi::new(bloom_size, map_length as u64, num_hash);
    let mut ref_kmer = fnv::FnvHashMap::default();
    let mut colors_accession = fnv::FnvHashMap::default();
    let mut accession_colors = fnv::FnvHashMap::default();
    let mut color = 0;
    let mut processed = 0;
    for accession in map.keys() {
        colors_accession.insert(color as usize, accession.to_owned()); //insert color, accession
        accession_colors.insert(accession, color); //insert accession, color
        color += 1; //color+= 1
    }
    let mut vec = Vec::with_capacity(batch);
    for (accession, v) in map {
        let color = accession_colors[accession];
        vec.push((color, v[0].to_owned()));
        if vec.len() % batch == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.1.to_owned());
                    let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m, 1);
                    (r.0.to_owned(), kmers)
                })
                .collect();
            //single threat
            for a in out_vec {
                let accession = a.0;
                let kmers = a.1;
                ref_kmer.insert(colors_accession[&accession].to_owned(), kmers.len());
                for k in &kmers {
                    bigsi.insert(accession as u64, k);
                }
            }
            processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    //one last time for the remainder of the accessions
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let fasta_vec = kmer::read_fasta(r.1.to_owned());
            let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m, 1);
            (r.0.to_owned(), kmers)
        })
        .collect();
    //single threat
    for a in out_vec {
        let accession = a.0;
        let kmers = a.1;
        ref_kmer.insert(colors_accession[&accession].to_owned(), kmers.len());
        for k in &kmers {
            bigsi.insert(accession as u64, k);
        } //
    }
    processed += vec.len();
    eprint!("processed {}/{} accessions\r", processed, map_length);
    bigsi.slim();
    (bigsi, colors_accession, ref_kmer)
}

pub fn build_multi_bigvec_mutex(
    map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    bloom_size: usize,
    num_hash: u64,
    k_size: usize,
    t: usize,
    quality: u8,
    cutoff: isize,
    batch: usize,
) -> (
    bigsi_rs::Bigsi,
    fnv::FnvHashMap<usize, String>,
    fnv::FnvHashMap<String, usize>,
) {
    //get a hash map with taxa and bit vector from bloom filter
    rayon::ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .unwrap();
    let map_length = map.len();
    //let mut bigsi = bi::Bigsi::new(bloom_size, map_length as u64, num_hash);
    let bigsi = Mutex::new(bigsi_rs::Bigsi::new(
        bloom_size,
        map_length as u64,
        num_hash,
    ));
    let ref_kmer = Mutex::new(fnv::FnvHashMap::default());
    let mut colors_accession = fnv::FnvHashMap::default();
    let mut accession_colors = fnv::FnvHashMap::default();
    let mut color = 0;
    let mut processed = 0;
    for accession in map.keys() {
        colors_accession.insert(color as usize, accession.to_owned()); //insert color, accession
        accession_colors.insert(accession, color); //insert accession, color
        color += 1; //color+= 1
    }
    let mut vec = Vec::with_capacity(batch);
    let mut map_vec = Vec::with_capacity(map_length);
    for (accession, v) in map {
        map_vec.push((accession, v));
    }
    for (accession, v) in map_vec {
        let color = accession_colors[accession];
        vec.push((color, v));
        if vec.len() % batch == 0 {
            processed += vec.len();
            //let mut out_vec: Vec<_> = vec![];
            //vec.par_iter()
            //    .map(|l| {
            vec.par_iter_mut().for_each(|l| {
                if l.1.len() == 2 {
                    let mut unfiltered =
                        kmer::kmers_fq_pe_qual(vec![&l.1[0], &l.1[1]], k_size, 1, quality);
                    let kmers = if cutoff == -1 {
                        let auto_cutoff = kmer::auto_cutoff(&unfiltered);
                        //let mut unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                        let keys_for_removal = kmer::get_removal_values(&unfiltered, auto_cutoff);
                        kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                        unfiltered
                    //kmer::clean_map(unfiltered, auto_cutoff)
                    } else {
                        let keys_for_removal =
                            kmer::get_removal_values(&unfiltered, cutoff as usize);
                        kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                        unfiltered
                        //kmer::clean_map(unfiltered, cutoff as usize)
                    };
                    let mut ref_kmer = ref_kmer.lock().unwrap();
                    ref_kmer.insert(colors_accession[&l.0].to_owned(), kmers.len());
                    let mut bigsi = bigsi.lock().unwrap();
                    for k in kmers.keys() {
                        bigsi.insert(l.0 as u64, k);
                    }
                //(l.0, kmers)
                } else {
                    if l.1[0].ends_with("gz") {
                        let mut unfiltered =
                            kmer::kmers_from_fq_qual(l.1[0].to_owned(), k_size, 1, 15);
                        let kmers = if cutoff == -1 {
                            let auto_cutoff = kmer::auto_cutoff(&unfiltered);
                            let keys_for_removal =
                                kmer::get_removal_values(&unfiltered, auto_cutoff);
                            kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                            unfiltered
                        //kmer::clean_map(unfiltered, auto_cutoff)
                        } else {
                            let keys_for_removal =
                                kmer::get_removal_values(&unfiltered, cutoff as usize);
                            kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                            unfiltered
                            //kmer::clean_map(unfiltered, cutoff as usize)
                        };
                        let mut ref_kmer = ref_kmer.lock().unwrap();
                        ref_kmer.insert(colors_accession[&l.0].to_owned(), kmers.len());
                        let mut bigsi = bigsi.lock().unwrap();
                        for k in kmers.keys() {
                            bigsi.insert(l.0 as u64, k);
                        }
                    //(l.0, kmers)
                    } else {
                        let vec = kmer::read_fasta(l.1[0].to_string());
                        let kmers = if cutoff == -1 {
                            kmer::kmerize_vector(&vec, k_size, 1)
                        } else {
                            let mut unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                            let keys_for_removal =
                                kmer::get_removal_values(&unfiltered, cutoff as usize);
                            kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                            unfiltered
                            //kmer::clean_map(unfiltered, cutoff as usize)
                        };
                        let mut ref_kmer = ref_kmer.lock().unwrap();
                        ref_kmer.insert(colors_accession[&l.0].to_owned(), kmers.len());
                        let mut bigsi = bigsi.lock().unwrap();
                        for k in kmers.keys() {
                            bigsi.insert(l.0 as u64, k);
                        }
                        //(l.0, kmers)
                    }
                }
            });
            //single threat
            /*
            for a in out_vec {
                let accession = a.0;
                let kmers = a.1;
                ref_kmer.insert(colors_accession[&accession].to_owned(), kmers.len());
                let mut bigsi = bigsi.lock().unwrap();
                for k in kmers.keys() {
                    bigsi.insert(accession as u64, k);
                }
            }*/
            //processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    //one last time for the remainder of the accessions
    processed += vec.len();
    vec.par_iter_mut().for_each(|l| {
        if l.1.len() == 2 {
            let mut unfiltered = kmer::kmers_fq_pe_qual(vec![&l.1[0], &l.1[1]], k_size, 1, quality);
            let kmers = if cutoff == -1 {
                let auto_cutoff = kmer::auto_cutoff(&unfiltered);
                //let mut unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                let keys_for_removal = kmer::get_removal_values(&unfiltered, auto_cutoff);
                kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                unfiltered
            //kmer::clean_map(unfiltered, auto_cutoff)
            } else {
                //let mut unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                let keys_for_removal = kmer::get_removal_values(&unfiltered, cutoff as usize);
                kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                unfiltered
                //kmer::clean_map(unfiltered, cutoff as usize)
            };
            let mut ref_kmer = ref_kmer.lock().unwrap();
            ref_kmer.insert(colors_accession[&l.0].to_owned(), kmers.len());
            let mut bigsi = bigsi.lock().unwrap();
            for k in kmers.keys() {
                bigsi.insert(l.0 as u64, k);
            }
        //(l.0, kmers)
        } else {
            if l.1[0].ends_with("gz") {
                let unfiltered = kmer::kmers_from_fq_qual(l.1[0].to_owned(), k_size, 1, 15);
                let kmers = if cutoff == -1 {
                    let auto_cutoff = kmer::auto_cutoff(&unfiltered);
                    kmer::clean_map(unfiltered, auto_cutoff)
                } else {
                    kmer::clean_map(unfiltered, cutoff as usize)
                };
                let mut ref_kmer = ref_kmer.lock().unwrap();
                ref_kmer.insert(colors_accession[&l.0].to_owned(), kmers.len());
                let mut bigsi = bigsi.lock().unwrap();
                for k in kmers.keys() {
                    bigsi.insert(l.0 as u64, k);
                }
            //(l.0, kmers)
            } else {
                let vec = kmer::read_fasta(l.1[0].to_string());
                let kmers = if cutoff == -1 {
                    kmer::kmerize_vector(&vec, k_size, 1)
                } else {
                    let mut unfiltered = kmer::kmerize_vector(&vec, k_size, 1);
                    let keys_for_removal = kmer::get_removal_values(&unfiltered, cutoff as usize);
                    kmer::clean_map_inplace(&mut unfiltered, &keys_for_removal);
                    unfiltered
                };
                let mut ref_kmer = ref_kmer.lock().unwrap();
                ref_kmer.insert(colors_accession[&l.0].to_owned(), kmers.len());
                let mut bigsi = bigsi.lock().unwrap();
                for k in kmers.keys() {
                    bigsi.insert(l.0 as u64, k);
                }
                //(l.0, kmers)
            }
        }
    });
    eprint!("processed {}/{} accessions\r", processed, map_length);
    let mut bigsi = bigsi.into_inner().unwrap();
    bigsi.slim();
    let ref_kmer = ref_kmer.into_inner().unwrap();
    (bigsi, colors_accession, ref_kmer)
}
