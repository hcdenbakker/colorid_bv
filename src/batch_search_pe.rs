use super::kmer;
use super::reports;
use bv::BitVec;
use bv::Bits;
use bv::BitsExt;
use fnv;
use std::collections::HashMap;
use std::time::SystemTime;

#[inline]
pub fn bitwise_and(vector_of_bitvectors: &Vec<&BitVec<u64>>) -> BitVec<u64> {
    let mut first = vector_of_bitvectors[0].to_owned();
    if vector_of_bitvectors.len() == 1 {
        first
    } else {
        for i in 1..vector_of_bitvectors.len() {
            let j = i as usize;
            first = first.bit_and(vector_of_bitvectors[j]).to_bit_vec();
        }
        first
    }
}

pub fn batch_search(
    files1: Vec<&str>,
    files2: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    n_ref_kmers: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
    k_size: usize,
    filter: isize,
    cov: f64,
    gene_search: bool,
    qual_offset: u8,
) {
    for (i, file1) in files1.iter().enumerate() {
        if file1.ends_with("gz") {
            let unfiltered = if files2.len() == 0 {
                eprintln!("{}", file1);
                eprintln!("Counting k-mers, this may take a while!");
                kmer::kmers_from_fq_qual(file1.to_owned().to_string(), k_size, 1, qual_offset)
            } else {
                eprintln!("Paired end: {} {}", file1, files2[i]);
                eprintln!("Counting k-mers, this may take a while!");
                kmer::kmers_fq_pe_qual(vec![&file1, &files2[i]], k_size, 1, qual_offset)
            };
            let kmers_query = if filter < 0 {
                let cutoff = kmer::auto_cutoff(&unfiltered);
                kmer::clean_map(unfiltered, cutoff)
            } else {
                kmer::clean_map(unfiltered, filter as usize)
            };
            let num_kmers = kmers_query.len() as f64;
            eprintln!("{} k-mers in query", num_kmers);
            let bigsi_search = SystemTime::now();
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for k in kmers_query.keys() {
                let bitvec = bigsi_map.get_bv(&k);
                /*let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index =
                        fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
                    if !bigsi_map.contains_key(&bit_index) {
                        break;
                    } else {
                        kmer_slices.push(&bigsi_map[&bit_index]);
                    }
                }*/
                if bitvec.is_empty() {
                    //if kmer_slices.len() < num_hash as usize {
                    continue;
                } else {
                    //let mut first = bitwise_and(&kmer_slices);
                    let mut hits = Vec::new();
                    for i in 0..bitvec.len() {
                        if bitvec[i] {
                            hits.push(&colors_accession[&(i as usize)]);
                        }
                    }
                    for h in &hits {
                        let count = report.entry(h.to_string()).or_insert(0);
                        *count += 1;
                    }
                    if hits.len() == 1 {
                        let key = hits[0];
                        let value = kmers_query[&k.to_string()] as f64;
                        uniq_freqs
                            .entry(key.to_string())
                            .or_insert_with(Vec::new)
                            .push(value);
                    }
                }
            }
            match bigsi_search.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Search: {} sec", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if !gene_search {
                reports::generate_report(
                    file1,
                    report,
                    &uniq_freqs,
                    &n_ref_kmers,
                    num_kmers as usize,
                    cov,
                );
            } else {
                reports::generate_report_gene(file1, report, num_kmers as usize, cov);
            }
        } else {
            //else we assume it is a fasta formatted file!
            eprintln!("{}", file1);
            eprintln!("Counting k-mers, this may take a while!");
            let vec_query = kmer::read_fasta(file1.to_owned().to_string());
            let unfiltered = kmer::kmerize_vector(&vec_query, k_size, 1);
            let kmers_query = if gene_search {
                kmer::clean_map(unfiltered, 0)
            } else if filter < 0 {
                eprintln!("no gene search");
                let cutoff = kmer::auto_cutoff(&unfiltered);
                kmer::clean_map(unfiltered, cutoff)
            } else {
                kmer::clean_map(unfiltered, filter as usize)
            };
            let num_kmers = kmers_query.len() as f64;
            eprintln!("{} k-mers in query", num_kmers);
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            let bigsi_search = SystemTime::now();
            for k in kmers_query.keys() {
                let bitvec = bigsi_map.get_bv(&k);
                /*let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index =
                        fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
                    if !bigsi_map.contains_key(&bit_index) {
                        break;
                    } else {
                        kmer_slices.push(&bigsi_map[&bit_index]);
                    }
                }*/
                if bitvec.is_empty() {
                    //if kmer_slices.len() < num_hash as usize {
                    continue;
                } else {
                    //let mut first = bitwise_and(&kmer_slices);
                    let mut hits = Vec::new();
                    for i in 0..bitvec.len() {
                        if bitvec[i] {
                            hits.push(&colors_accession[&(i as usize)]);
                        }
                    }
                    for h in &hits {
                        let count = report.entry(h.to_string()).or_insert(0);
                        *count += 1;
                    }
                    if hits.len() == 1 {
                        let key = hits[0];
                        let value = kmers_query[&k.to_string()] as f64;
                        uniq_freqs
                            .entry(key.to_string())
                            .or_insert_with(Vec::new)
                            .push(value);
                    }
                }
            }
            match bigsi_search.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Search: {} sec", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if !gene_search {
                reports::generate_report(
                    file1,
                    report,
                    &uniq_freqs,
                    &n_ref_kmers,
                    num_kmers as usize,
                    cov,
                );
            } else {
                reports::generate_report_gene(file1, report, num_kmers as usize, cov);
            }
        }
    }
}

pub fn batch_search_mf(
    files: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    n_ref_kmers: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
    k_size: usize,
    cov: f64,
) {
    for file in files {
        let (labels, sequences) = super::kmer::read_fasta_mf(file.to_owned());
        for (i, label) in labels.iter().enumerate() {
            let kmers_query = super::kmer::kmerize_string(sequences[i].to_owned(), k_size);
            match kmers_query{
            None    => println!("Warning! no kmers in query '{}'; maybe your kmer length is larger than your query length?", label),
            Some(kmers_query) =>
            {
            let num_kmers = kmers_query.len() as f64;
            let mut report = HashMap::new();
            let mut uniq_freqs = HashMap::new();
            for k in kmers_query.keys() {
                let bitvec = bigsi_map.get_bv(&k);
                /*let mut kmer_slices = Vec::new();
                for i in 0..num_hash {
                    let bit_index =
                        fasthash::xx::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
                    if !bigsi_map.contains_key(&bit_index) {
                        break;
                    } else {
                        kmer_slices.push(&bigsi_map[&bit_index]);
                    }
                }*/
                if bitvec.is_empty() {
                    //if kmer_slices.len() < num_hash as usize {
                    continue;
                } else {
                    //let mut first = bitwise_and(&kmer_slices);
                    let mut hits = Vec::new();
                    for i in 0..bitvec.len() {
                        if bitvec[i] {
                            hits.push(&colors_accession[&(i as usize)]);
                        }
                    }
                    for h in &hits {
                        let count = report.entry(h.to_string()).or_insert(0);
                        *count += 1;
                    }
                    if hits.len() == 1 {
                        let key = hits[0];
                        let value = kmers_query[&k.to_string()] as f64;
                        uniq_freqs
                            .entry(key.to_string())
                            .or_insert_with(Vec::new)
                            .push(value);
                    }
                }
            }
            reports::generate_report_gene(label, report, num_kmers as usize, cov);
        },
        }
        }
    }
}
