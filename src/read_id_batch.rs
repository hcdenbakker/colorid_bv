use super::build;
use super::read_id_mt_pe;
use super::reports;
use bincode::deserialize_from;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

#[allow(unused_assignments)]
pub fn read_id_batch(
    batch_samples: &str,
    index: &str,
    down_sample: usize,
    fp_correct: f64,
    batch: usize,
    quality: u8,
    bitvector_sample: usize,
    high_mem_load: bool,
    tag: &str,
    supress_taxon: &str,
) {
    //tab to map to read batch
    let batch_map = build::tab_to_map(batch_samples.to_string());
    let now = Instant::now();
    eprintln!("Loading index...");
    let mut db = bigsi_rs::Bigsi::default();
    db.read(&(index.to_owned() + "/bigsi"));

    //parameters
    let mut reader = BufReader::new(
        File::open(&(index.to_owned() + "/parameters")).expect("Can't open parameters!"),
    );
    let parameters: build::Parameters = deserialize_from(&mut reader).expect("can't deserialize");
    //colors
    let mut reader = BufReader::new(
        File::open(&(index.to_owned() + "/colors")).expect("Can't open parameters!"),
    );
    let colors = deserialize_from(&mut reader).expect("can't deserialize");
    //ref_kmers
    let mut reader = BufReader::new(
        File::open(&(index.to_owned() + "/ref_kmers")).expect("Can't open parameters!"),
    );
    let ref_kmers = deserialize_from(&mut reader).expect("can't deserialize");
    eprintln!("{}", now.elapsed().as_secs());

    //start iterating through files here
    for (accession, fq_vec) in &batch_map {
        eprintln!("Classifying {}", accession);
        let prefix = format!("{}_{}", accession, tag);
        let mut fq: Vec<_> = vec![];
        fq = fq_vec.iter().map(|r| &r[..]).collect();
        if fq[0].ends_with(".gz") {
            if fq.len() > 1 {
                read_id_mt_pe::per_read_stream_pe_itertools(
                    fq,
                    &db,
                    &colors,
                    &ref_kmers,
                    db.bigsi.len() as u64,
                    db.num_hashes,
                    parameters.k_size,
                    parameters.m_size,
                    down_sample,
                    fp_correct,
                    batch,
                    &prefix,
                    quality,
                    bitvector_sample,
                    supress_taxon,
                )
            } else {
                read_id_mt_pe::per_read_stream_se(
                    fq,
                    &db,
                    &colors,
                    &ref_kmers,
                    db.bigsi.len() as u64,
                    db.num_hashes,
                    parameters.k_size,
                    parameters.m_size,
                    down_sample,
                    fp_correct,
                    batch,
                    &prefix,
                    quality,
                    bitvector_sample,
                    supress_taxon,
                )
            };
        } else {
            read_id_mt_pe::stream_fasta(
                fq,
                &db,
                &colors,
                &ref_kmers,
                db.bigsi.len() as u64,
                db.num_hashes,
                parameters.k_size,
                parameters.m_size,
                down_sample,
                fp_correct,
                batch,
                &prefix,
                bitvector_sample,
                supress_taxon,
            );
        }
        reports::read_counts_five_fields(prefix.to_owned() + "_reads.txt", &prefix);
    }
}
