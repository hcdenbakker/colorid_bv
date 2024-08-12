use super::kmer;
use super::seq;
use bv::BitVec;
use bv::*;
use xxh3;
use flate2::read::MultiGzDecoder;
use fnv;
use itertools::Itertools;
use probability::prelude::*;
use rayon::prelude::*;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::sync::Arc;
use std::time::SystemTime;

pub fn false_prob_map(
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
) -> fnv::FnvHashMap<usize, f64> {
    //first change a colors to accessions map to a accessions to colors map
    let mut accession_color = HashMap::new();
    for (color, accession) in colors_accession {
        accession_color.insert(accession, color);
    }
    let mut false_positive_p = fnv::FnvHashMap::default();
    for (key, value) in ref_kmers_in {
        let color = accession_color.get(key).unwrap();
        false_positive_p.insert(
            *color.to_owned(),
            false_prob(bloom_size as f64, num_hash as f64, *value as f64),
        );
    }
    false_positive_p
}

#[inline]
pub fn bitwise_and(vector_of_bitvectors: &[&BitVec<u64>]) -> BitVec<u64> {
    let mut first = vector_of_bitvectors[0].to_owned();
    if vector_of_bitvectors.len() == 1 {
        first
    } else {
        for i in 1..vector_of_bitvectors.len() {
            let j = i as usize;
            first = first.bit_and(&vector_of_bitvectors[j]).to_bit_vec();
        }
        first
    }
}

/*
#[inline]
pub fn bitwise_and(vector_of_bitvectors: &Vec<&BitVec>) -> BitVec {
    if vector_of_bitvectors.len() == 1 {
        vector_of_bitvectors[0].to_owned()
    } else {
        let first = vector_of_bitvectors[0].bit_and(&vector_of_bitvectors[1]);
        if vector_of_bitvectors.len() > 2 {
            for i in 2..vector_of_bitvectors.len() {
                let j = i as usize;
                let first = first.bit_and(&vector_of_bitvectors[j]);
            }
        }
        first.to_bit_vec()
    }
}*/

//change a vector of strings to a comma delimited string
#[inline]
pub fn vec_strings_to_string(vector_in: &[String]) -> String {
    let mut comma_separated = String::new();
    for s in vector_in {
        comma_separated.push_str(&s.to_string());
        comma_separated.push_str(",");
    }
    comma_separated.pop();
    comma_separated
}

pub fn search_index_classic(
    bigsi_map: &fnv::FnvHashMap<u64, BitVec<u64>>,
    map: &fnv::FnvHashSet<std::string::String>,
    bloom_size: u64,
    num_hash: u64,
    no_hits_num: usize,
) -> fnv::FnvHashMap<usize, usize> {
    let mut report = fnv::FnvHashMap::default();
    let empty_bitvec: BitVec<u64> = BitVec::new_fill(false, no_hits_num as u64);
    for k in map {
        let mut kmer_slices = Vec::new();
        for i in 0..num_hash {
            let bit_index =
                xxh3::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
            //let bi = bit_index as usize;
            if !bigsi_map.contains_key(&bit_index) || bigsi_map[&bit_index] == empty_bitvec {
                break;
            } else {
                kmer_slices.push(bigsi_map.get(&bit_index).unwrap());
            }
        }
        if kmer_slices.len() < num_hash as usize {
            *report.entry(no_hits_num).or_insert(0) += 1;
            break;
        } else {
            let first = bitwise_and(&kmer_slices);
            //let mut color = 0;
            for item in 0..first.len() {
                if first[item] {
                    *report.entry(item as usize).or_insert(0) += 1;
                }
                //  color += 1;
            }
        }
    }
    report
}

#[inline]
pub fn search_index_bigvec(
    bigsi_map: &bigsi_rs::Bigsi,
    map: &fnv::FnvHashSet<std::string::String>,
    no_hits_num: usize,
    start_sample: usize,
) -> fnv::FnvHashMap<usize, usize> {
    let mut report = fnv::FnvHashSet::default();
    let mut final_report = fnv::FnvHashMap::default();
    let mut counter = 0;
    for k in map {
        if start_sample == 0 || counter < start_sample {
            let first = bigsi_map.get_bv(&k);
            if first.is_empty() {
                *final_report.entry(no_hits_num).or_insert(0) += 1;
            } else {
                //let first = bitwise_and(&kmer_slices);
                //let mut color: usize = 0;
                for item in 0..first.len() {
                    if first[item] {
                        report.insert(item);
                        *final_report.entry(item as usize).or_insert(0) += 1;
                    }
                    //color += 1;
                }
            }
        } else {
            let first = bigsi_map.get_bv(&k);
            if first.is_empty() {
                *final_report.entry(no_hits_num).or_insert(0) += 1;
                break;
            } else {
                for item in &report {
                    if first[*item] {
                        *final_report.entry(*item as usize).or_insert(0) += 1;
                    }
                }
            }
        }
        counter += 1;
    }
    final_report
}

#[inline]
pub fn not_fp_signicant(
    observations: usize,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    fp_correct: f64,
    taxon: usize,
    taxon_hits: usize,
) -> bool {
    let p_false = child_fp.get(&taxon).unwrap();
    let distribution = Binomial::new(observations, *p_false);
    let critical_value = observations as f64 * p_false;
    let mpf = distribution.mass(taxon_hits);
    ((taxon_hits as f64) < critical_value)
        || (((taxon_hits as f64) > critical_value) && (mpf >= fp_correct))
}

//kmer poll classification plus raw count data as output
// 1. filter hits below false prob threshold
// 2. find tophits
// 3 if only one tophit, report with unique flag
pub fn kmer_poll_plus<'a>(
    report: &fnv::FnvHashMap<usize, usize>,
    kmer_length: usize,
    colors_accession: &'a fnv::FnvHashMap<usize, String>,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    no_hits_num: usize,
    fp_correct: f64,
    supress_accession: usize,
) -> (String, usize, usize, &'a str, usize) {
    let mut count_vec: Vec<_> = report.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    if (count_vec[0].0 == &no_hits_num) && (count_vec.len() == 1) {
        return (
            "no_hits".to_string(),
            0 as usize,
            kmer_length,
            "accept",
            0 as usize,
        );
    };
    let mut significant_hits = Vec::new();
    for t in &count_vec {
        if t.0 == &no_hits_num
            || t.0 == &supress_accession
            || not_fp_signicant(kmer_length, child_fp, fp_correct, *t.0, *t.1)
        {
            continue;
        } else {
            significant_hits.push(t.to_owned());
        }
    }
    if significant_hits.is_empty() {
        (
            "no_significant_hits".to_string(),
            0 as usize,
            kmer_length,
            "reject",
            0 as usize,
        )
    } else {
        //significant_hits.sort_by(|a, b| b.1.cmp(a.1));
        let first_tophit = colors_accession[&significant_hits[0].0].to_owned();
        let mut top_hits = Vec::new();
        for h in &significant_hits {
            if h.1 == significant_hits[0].1 {
                top_hits.push(colors_accession[&h.0].to_owned())
            }
        }
        if top_hits.len() == 1 {
            (
                first_tophit,
                significant_hits[0].1.to_owned(),
                kmer_length,
                "accept",
                top_hits.len(),
            )
        } else {
            (
                vec_strings_to_string(&top_hits),
                significant_hits[0].1.to_owned(),
                kmer_length,
                "reject",
                top_hits.len(),
            )
        }
    }
}

#[derive(Default)]
pub struct SeqRead {
    pub id: String,       //id including >
    pub seq: Vec<String>, // sequence
}

impl SeqRead {
    pub fn new() -> SeqRead {
        SeqRead {
            id: String::new(),
            seq: Vec::new(),
        }
    }
}

#[derive(Default)]
pub struct SeqReadstr<'a> {
    pub id: &'a str,       //id including >
    pub seq: Vec<&'a str>, // sequence
}

impl<'a> SeqReadstr<'a> {
    pub fn new() -> SeqReadstr<'a> {
        SeqReadstr {
            id: "",
            seq: Vec::new(),
        }
    }
}

#[allow(unused_assignments)]
pub fn parallel_vec_bigvec(
    vec: &[(String, Vec<String>)],
    bigsi: &bigsi_rs::Bigsi,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
    no_hits_num: usize,
    k: usize,
    m: usize,
    d: usize,
    fp_correct: f64,
    start_sample: usize,
    supress_taxon: usize,
) -> std::vec::Vec<(std::string::String, String, usize, usize, String, usize)> {
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&bigsi_rs::Bigsi> = Arc::new(bigsi);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if r.1[0].len() < k {
                (
                    r.0.to_owned(),
                    "too_short".to_string(),
                    0 as usize,
                    0 as usize,
                    "accept".to_string(),
                    0 as usize,
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n_set(&r.1, k, d)
                } else {
                    kmer::minimerize_vector_skip_n_set(&r.1, k, m, d)
                };
                let report = //if start_sample == 0 {
                    //search_index_classic(&child_bigsi, &map, bloom_size, num_hash, no_hits_num)
                //} else {
                    search_index_bigvec(
                        &child_bigsi,
                        &map,
                        no_hits_num,
                        start_sample,
                    )
                ;
                if report.is_empty() {
                    (
                        r.0.to_owned(),
                        "no_hits".to_string(),
                        0 as usize,
                        map.len(),
                        "accept".to_string(),
                        0 as usize,
                    )
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        map.len(),
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                        supress_taxon,
                    );
                    (
                        r.0.to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3.to_owned(),
                        classification.4,
                    )
                }
            }
        })
        .collect();
    out_vec
}
/*
#[allow(unused_assignments)]
pub fn parallel_vec_str(
    vec: &Vec<(&str, Vec<&str>)>,
    bigsi: &fnv::FnvHashMap<usize, BitVec>,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
    k: usize,
    m: usize,
    d: usize,
    fp_correct: f64,
    start_sample: usize,
) -> std::vec::Vec<(std::string::String, String, usize, usize, String, usize)> {
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&fnv::FnvHashMap<usize, BitVec>> = Arc::new(bigsi);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if r.1[0].len() < k {
                (
                    r.0.to_owned(),
                    "too_short".to_string(),
                    0 as usize,
                    0 as usize,
                    "accept".to_string(),
                    0 as usize,
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n_set_str(&r.1, k, d)
                } else {
                    kmer::minimerize_vector_skip_n_set_str(&r.1, k, m, d)
                };
                let report = if start_sample == 0 {
                    search_index_classic(&child_bigsi, &map, bloom_size, num_hash, no_hits_num)
                } else {
                    search_index(
                        &child_bigsi,
                        &map,
                        bloom_size,
                        num_hash,
                        no_hits_num,
                        start_sample,
                    )
                };
                if report.is_empty() {
                    (
                        r.0.to_owned(),
                        "no_hits".to_string(),
                        0 as usize,
                        map.len(),
                        "accept".to_string(),
                        0 as usize,
                    )
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        map.len(),
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    );
                    (
                        r.0.to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3.to_owned(),
                        classification.4,
                    )
                }
            }
        })
        .collect();
    out_vec
}*/

#[allow(unused_assignments)]
pub fn stream_fasta(
    filenames: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
    k: usize,
    m: usize,
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    start_sample: usize,
    supress_accession: &str,
) {
    let f = File::open(filenames[0]).expect("file not found");
    //let iter1 = io::BufReader::new(f).lines();
    let mut iter1 = io::BufReader::new(f);
    let mut l = String::with_capacity(250);
    let mut vec = Vec::with_capacity(b);
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let mut sub_string = String::new();
    let search_time = SystemTime::now();
    let mut count = 0;
    let mut read_count = 0;
    let mut fasta = SeqRead::new();
    let supress_taxon: usize = if supress_accession == "None" {
        colors_accession.len() + 1
    } else {
        //we need to lookup the color of the taxon to be surpressed
        let mut surpress = colors_accession.len() + 1;
        let target = supress_accession.to_string();
        for (k, v) in colors_accession {
            if v == &target {
                surpress = *k;
                break;
            }
        }
        surpress
    };
    //for line in iter1 {
    while iter1.read_line(&mut l).unwrap() > 0 {
        //let l = line.unwrap();
        if count == 0 {
            fasta.id = l[..l.len() - 1].to_string();
        } else {
            if l.contains('>') {
                if !sub_string.is_empty() {
                    fasta.seq = vec![sub_string.to_string()];
                    vec.push((fasta.id, fasta.seq));
                    fasta.id = l[..l.len() - 1].to_string();
                    sub_string.clear();
                }
            } else {
                sub_string.push_str(&l);
            }
        }
        count += 1;
        if vec.len() % b == 0 {
            let c = parallel_vec_bigvec(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
                supress_taxon,
            );
            read_count += c.len();
            eprint!(" {} reads classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
        l.clear();
    }
    fasta.seq = vec![sub_string];
    vec.push((fasta.id, fasta.seq));
    let c = parallel_vec_bigvec(
        &vec,
        bigsi_map,
        colors_accession,
        ref_kmers_in,
        bloom_size,
        num_hash,
        no_hits_num,
        k,
        m,
        d,
        fp_correct,
        start_sample,
        supress_taxon,
    );
    read_count += c.len();
    for id in c {
        file.write_all(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                id.0, id.1, id.2, id.3, id.4, id.5
            )
            .as_bytes(),
        )
        .expect("could not write results!");
    }
    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}
/*
#[allow(unused_assignments)]
pub fn stream_fasta_str(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    start_sample: usize,
)
{
    let f = File::open(filenames[0]).expect("file not found");
    let iter1 = io::BufReader::new(f).lines();
    let mut vec = Vec::with_capacity(b);
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let mut sub_string = String::new();
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .expect("Can't initialize ThreadPoolBuilder");
    let search_time = SystemTime::now();
    let mut count = 0;
    let mut read_count = 0;
    let mut fasta = SeqReadstr::new();
    for line in iter1 {
        let l = line.unwrap();
        if count == 0 {
            fasta.id = &l;
        } else {
            if l.contains('>') {
                if sub_string.len() > 0 {
                    fasta.seq = vec![&sub_string];
                    vec.push((fasta.id, fasta.seq));
                    fasta.id = &l;
                    sub_string.clear();
                }
            } else {
                sub_string.push_str(&l);
            }
        }
        count += 1;
        if vec.len() % b == 0 {
            let c = parallel_vec_str(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
            );
            read_count += c.len();
            eprint!(" {} reads classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
    }
    fasta.seq = vec![&sub_string];
    vec.push((fasta.id, fasta.seq));
    let c = parallel_vec_str(
        &vec,
        bigsi_map,
        colors_accession,
        ref_kmers_in,
        bloom_size,
        num_hash,
        no_hits_num,
        k,
        m,
        d,
        fp_correct,
        start_sample,
    );
    read_count += c.len();
    for id in c {
        file.write_all(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                id.0, id.1, id.2, id.3, id.4, id.5
            )
            .as_bytes(),
        )
        .expect("could not write results!");
    }
    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}
*/

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k)
}

#[allow(unused_assignments)]
pub fn per_read_stream_pe_itertools(
    filenames: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi, //has to be an Arc
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
    k: usize,
    m: usize,
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    qual_offset: u8,
    start_sample: usize,
    supress_accession: &str,
) {
    let no_hits_num: usize = colors_accession.len();
    let supress_taxon: usize = if supress_accession == "None" {
        colors_accession.len() + 1
    } else {
        //we need to lookup the color of the taxon to be surpressed
        let mut surpress = colors_accession.len() + 1;
        let target = supress_accession.to_string();
        for (k, v) in colors_accession {
            if v == &target {
                surpress = *k;
                break;
            }
        }
        surpress
    };
    let search_time = SystemTime::now();
    let mut read_count = 0;
    let mut vec = Vec::with_capacity(b);
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let f1 = File::open(&filenames[0]).expect("file not found");
    let f2 = File::open(&filenames[1]).expect("file not found");
    let d1 = MultiGzDecoder::new(f1);
    let d2 = MultiGzDecoder::new(f2);
    let iter1 = BufReader::new(d1).lines().chunks(b * 4);
    let iter2 = BufReader::new(d2).lines().chunks(b * 4);
    for (chunk1, chunk2) in izip!(&iter1, &iter2) {
        for ((header, sequence, _, qual), (_, sequence2, _, qual2)) in
            izip!(chunk1.tuples(), chunk2.tuples())
        {
            //total_bp += (sequence.unwrap().len() + sequence2.unwrap().len());
            //count += 1;
            vec.push((
                header.unwrap(),
                vec![
                    super::seq::qual_mask(&sequence.unwrap(), &qual.unwrap(), qual_offset),
                    super::seq::qual_mask(&sequence2.unwrap(), &qual2.unwrap(), qual_offset),
                ],
            ));
        }
        let c = parallel_vec_bigvec(
            &vec,
            bigsi_map,
            colors_accession,
            ref_kmers_in,
            bloom_size,
            num_hash,
            no_hits_num,
            k,
            m,
            d,
            fp_correct,
            start_sample,
            supress_taxon,
        );
        read_count += c.len();
        eprint!("{} read pairs classified\r", read_count);
        for id in c {
            file.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    id.0, id.1, id.2, id.3, id.4, id.5
                )
                .as_bytes(),
            )
            .expect("could not write results!");
        }
        vec.clear();
    }
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} read pairs in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}

#[allow(unused_assignments)]
pub fn per_read_stream_se(
    filenames: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: u64,
    num_hash: u64,
    k: usize,
    m: usize, //0 == no m, otherwise minimizer
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    qual_offset: u8,
    start_sample: usize,
    supress_accession: &str,
) /* -> std::collections::HashMap<std::string::String, usize>*/
{
    let search_time = SystemTime::now();
    let mut vec = Vec::with_capacity(b);
    let mut line_count = 1;
    let _header = "".to_string();
    let _sequence = "".to_string();
    let no_hits_num: usize = colors_accession.len();
    let supress_taxon: usize = if supress_accession == "None" {
        colors_accession.len() + 1
    } else {
        //we need to lookup the color of the taxon to be surpressed
        let mut surpress = colors_accession.len() + 1;
        let target = supress_accession.to_string();
        for (k, v) in colors_accession {
            if v == &target {
                surpress = *k;
                break;
            }
        }
        surpress
    };
    let mut read_count = 0;
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let batch = b * 4;
    let _num_files = filenames.len();
    let f = File::open(&filenames[0]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let iter1 = io::BufReader::new(d1).lines();
    let mut fastq = seq::Fastq::new();
    for line in iter1 {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            fastq.id = l.to_string();
        } else if line_count % 4 == 2 {
            fastq.seq1 = l.to_string();
        } else if line_count % 4 == 0 {
            vec.push((
                fastq.id.to_owned(),
                vec![seq::qual_mask(&fastq.seq1, &l, qual_offset)],
            ));
        }
        line_count += 1;
        if line_count % batch == 0 {
            let c = parallel_vec_bigvec(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
                supress_taxon,
            );
            read_count += c.len();
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
    }
    //and the last block if there is any left
    let c = parallel_vec_bigvec(
        &vec,
        bigsi_map,
        colors_accession,
        ref_kmers_in,
        bloom_size,
        num_hash,
        no_hits_num,
        k,
        m,
        d,
        fp_correct,
        start_sample,
        supress_taxon,
    );

    read_count += c.len();
    eprint!("{} read pairs classified\r", read_count);
    for id in c {
        file.write_all(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                id.0, id.1, id.2, id.3, id.4, id.5
            )
            .as_bytes(),
        )
        .expect("could not write results!");
    }
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}
