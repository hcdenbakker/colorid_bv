use super::kmer;
use bv::BitVec;
use bv::Bits;
use bv::BitsExt;
use fnv;

#[inline]
pub fn bitwise_and(vector_of_bitvectors: &Vec<BitVec<usize>>) -> BitVec<usize> {
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

pub fn batch_search(
    files: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    k_size: usize,
    _cov: f64,
) {
    for file in files {
        //only fasta formatted file!
        eprintln!("Counting k-mers, this may take a while!");
        let vec_query = kmer::read_fasta(file.to_owned());
        let kmers_query = kmer::kmerize_vector(&vec_query, k_size, 1);
        eprintln!("{} kmers in query", kmers_query.len());
        if kmers_query.len() == 0 {
            eprintln!("Warning! no kmers in query; maybe your kmer length is larger than your query length?");
        } else {
            let mut kmer_slices = Vec::new();
            for k in kmers_query.keys() {
                let bitvec = bigsi_map.get_bv(&k);
                if bitvec.is_empty() {
                    break;
                } else {
                    kmer_slices.push(bitvec);
                }
            }
            if kmer_slices.len() < kmers_query.len() {
                eprintln!("No perfect hits!");
            } else {
                //bit-wise AND
                let first = bitwise_and(&kmer_slices);
                let mut hits = Vec::new();
                for i in 0..first.len() {
                    if first[i] {
                        hits.push(&colors_accession[&(i as usize)]);
                    }
                }
                eprintln!("{} hits", hits.len());
                for h in &hits {
                    println!("{}\t{}\t{}\t1.00", file, h, kmers_query.len());
                }
            }
        }
    }
}

pub fn batch_search_mf(
    files: Vec<&str>,
    bigsi_map: &bigsi_rs::Bigsi,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    k_size: usize,
    _cov: f64,
) {
    for file in files {
        let (labels, sequences) = super::kmer::read_fasta_mf(file.to_owned());
        for (i, label) in labels.iter().enumerate() {
            //only fasta formatted file!
            //eprintln!("Counting k-mers, this may take a while!");
            let kmers = super::kmer::kmerize_string(sequences[i].to_owned(), k_size);
            //let kmers_query = kmer_fa::clean_map(unfiltered, filter);
            match kmers{
                None    => println!("Warning! no kmers in query '{}'; maybe your kmer length is larger than your query length?", label),
                Some(kmers_query) =>
                {eprintln!("{} kmers in query", kmers_query.len());
                let mut kmer_slices = Vec::new();
                for k in kmers_query.keys() {
                let bitvec = bigsi_map.get_bv(&k);
                if bitvec.is_empty(){
                    break;
                    } else {
                        kmer_slices.push(bitvec);
                    }
                }
                if kmer_slices.len() < kmers_query.len() {
                    eprintln!("No perfect hits!");
                } else {
                    //bit-wise AND
                    let first = bitwise_and(&kmer_slices);
                    let mut hits = Vec::new();
                    for i in 0..first.len() {
                        if first[i] {
                            hits.push(&colors_accession[&(i as usize)]);
                        }
                    }
                    eprintln!("{} hits", hits.len());
                    for h in &hits {
                        println!("{}\t{}\t{}\t1.00", label.to_string(), h, kmers_query.len());
                    }
                }
                },
            }
        }
    }
}
