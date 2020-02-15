extern crate colorid_bigvec as colorid;
extern crate rayon;
#[macro_use]
extern crate clap;

use bincode::{deserialize_from, serialize};
use clap::{App, AppSettings, Arg, SubCommand};
//use colorid_bigvec::bi;
use colorid_bigvec::build;
use colorid_bigvec::read_id_mt_pe;
use colorid_bigvec::read_id_mt_pe::false_prob;
use rayon::ThreadPoolBuilder;
use std::alloc::System;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::time::SystemTime;

#[global_allocator]
static GLOBAL: System = System;

fn main() -> std::io::Result<()> {
    let matches = App::new("colorid")
        .version("0.1.4.2")
        .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
        .about("BIGSI based taxonomic ID of sequence data")
        .setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("build")
                .about("builds a bigsi")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true),
                )
                /*.help(
                              "                              -b, --bigsi=[FILE] 'Sets the prefix of the index file'
                              -r, --refs      'two column tab delimited file, first column taxon name, second column file  with sequence data plus path'
                              -k, --kmer      'sets kmer size to use for index'
                              -n, --num_hashes  'number of hashes to use for bloom filters'
                              -s, --bloom     'size bloom filter'
                              -m, --minimizer 'build index from minimizers'
                              -t, --threads 'sets number of threads, if set to 0, takes as many threads as it can get, default 1'")*/
                .arg(
                    Arg::with_name("ref_file")
                        .help("Sets the input file to use")
                        .required(true)
                        .short("r")
                        .takes_value(true)
                        .long("refs"),
                )
                .arg(
                    Arg::with_name("k-mer_size")
                        .help("Sets k-mer size")
                        .required(true)
                        .short("k")
                        .takes_value(true)
                        .long("kmer"),
                )
                .arg(
                    Arg::with_name("num_hashes")
                        .help("Sets number of hashes for bloom filter")
                        .required(true)
                        .short("n")
                        .takes_value(true)
                        .long("num_hashes"),
                )
                .arg(
                    Arg::with_name("length_bloom")
                        .help("Sets the length of the bloom filter")
                        .required(true)
                        .short("s")
                        .takes_value(true)
                        .long("bloom"),
                )
                .arg(
                    Arg::with_name("minimizer")
                        .help("build index with minimizers of set value")
                        .required(false)
                        .short("m")
                        .takes_value(true)
                        .long("minimizer"),
                )
                .arg(
                    Arg::with_name("threads")
                        .help("number of threads to use, if not set one thread will be used")
                        .required(false)
                        .short("t")
                        .takes_value(true)
                        .long("threads"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("minimum phred score to keep basepairs within read (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                )
                .arg(
                    Arg::with_name("batch")
                        .help("Sets size of batch of accessions to be processed in parallel (default 300)")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("batch"),
                )
                .arg(
                    Arg::with_name("filter")
                        .help("minimum coverage kmer threshold (default automatically detected)")
                        .required(false)
                        .short("f")
                        .takes_value(true)
                        .long("filter"),
                        ),
        )
        .subcommand(
            SubCommand::with_name("search")
                .about("does a bigsi search on one or more fasta/fastq.gz files")
                .version("0.1")
                .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .help("Sets the name of the index file for search")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true),
                )
                /*.help(
                              "                              -b, --bigsi=[FILE] 'Sets the name of the index file for search'
                              -q, --query      'one or more fastq.gz or fasta formatted files to be queried'
                              -r, --reverse    'one or more fastq.gz (not fasta!) formatted files, supposed to be reverse reads of a paired end run'  
                              -f, --filter     'Sets threshold to filter k-mers by frequency'
                              -p, --p_shared        'minimum proportion of kmers shared with reference'
                              -g, If set('-g'), the proportion of kmers from the query matching the entries in the index will be reported'
                              -s, If set('-s'), the 'speedy' 'perfect match' alforithm will be performed'
                              -m, If set('-m'), each accession in a multifasta will betreated as a separate query, currently only with the -s option'")*/
                .arg(
                    Arg::with_name("query")
                        .help("query file(-s)fastq.gz")
                        .required(true)
                        .min_values(1)
                        .short("q")
                        .takes_value(true)
                        .long("query"),
                )
                .arg(
                    Arg::with_name("reverse")
                        .help("reverse file(-s)fastq.gz")
                        .required(false)
                        .min_values(1)
                        .default_value("none")
                        .short("r")
                        .takes_value(true)
                        .long("reverse"),
                )
                .arg(
                    Arg::with_name("filter")
                        .help("set minimum k-mer frequency ")
                        .required(false)
                        .short("f")
                        .takes_value(true)
                        .long("filter"),
                )
                .arg(
                    Arg::with_name("shared_kmers")
                        .help("set minimum proportion of shared k-mers with a reference")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("p_shared"),
                )
                .arg(
                    Arg::with_name("gene_search")
                        .help("If set('-g'), the proportion of kmers from the query matching the entries in the index will be reported")
                        .required(false)
                        .short("g")
                        .takes_value(false)
                        .long("gene_search"),
                )
                .arg(
                    Arg::with_name("perfect_search")
                        .help("If ('-s') is set, the fast 'perfect match' algorithm will be used")
                        .required(false)
                        .short("s")
                        .takes_value(false)
                        .long("perfect_search"),
                )
                .arg(
                    Arg::with_name("multi_fasta")
                        .help("If set('-m'), each accession in a multifasta will betreated as a separate query, currently only with the -s option")
                        .required(false)
                        .short("m")
                        .takes_value(false)
                        .long("multi_fasta"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("minimum phred score to keep basepairs within read (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                ),
        )
        .subcommand(
            SubCommand::with_name("info")
                .about("dumps index parameters and accessions")
                .version("0.1")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true)
                        .help("index for which info is requested"),
                        )
                .arg(
                    Arg::with_name("compressed")
                        .help("If set to 'true', it is assumed a gz compressed index is used")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("compressed"),
                ),
        )
        .subcommand(
            SubCommand::with_name("read_id")
                .about("id's reads")
                .version("0.2")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true)
                        .help("index to be used for search"),
                        )
                .arg(
                    Arg::with_name("query")
                        .help("query file(-s)fastq.gz")
                        .required(true)
                        .min_values(1)
                        .short("q")
                        .takes_value(true)
                        .long("query"),
                )
                .arg(
                    Arg::with_name("batch")
                        .help("Sets size of batch of reads to be processed in parallel (default 50,000)")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("batch"),
                )
                .arg(
                    Arg::with_name("threads")
                        .help("number of threads to use, if not set the maximum available number threads will be used")
                        .required(false)
                        .short("t")
                        .takes_value(true)
                        .long("threads"),
                )
                .arg(
                    Arg::with_name("prefix")
                        .help("prefix for output file(-s)")
                        .required(true)
                        .short("n") //running out of options here!
                        .takes_value(true)
                        .long("prefix"),
                )
                .arg(
                    Arg::with_name("down_sample")
                        .help("down-sample k-mers used for read classification, default 1; increases speed at cost of decreased sensitivity ")
                        .required(false)
                        .short("d")
                        .takes_value(true)
                        .long("down_sample"),
                )
                .arg(
                    Arg::with_name("high_mem_load")
                        .help("When this flag is set, a faster, but less memory efficient method to load the index is used. Loading the index requires approximately 2X the size of the index of RAM. ")
                        .required(false)
                        .short("H")
                        .takes_value(false)
                        .long("high_mem_load"),
                )
                .arg(
                    Arg::with_name("fp_correct")
                        .help("Parameter to correct for false positives, default 3 (= 0.001), maybe increased for larger searches. Adjust for larger datasets")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("fp_correct"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("kmers with nucleotides below this minimum phred score will be excluded from the analyses (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                )
                .arg(
                    Arg::with_name("bitvector_sample")
                        .help("Collects matches for subset of kmers indicated (default=3), using this subset to more rapidly find hits for the remainder of the kmers")
                        .required(false)
                        .short("B")
                        .takes_value(true)
                        .long("bitvector_sample"),
                        )
                .arg(
                    Arg::with_name("supress_taxon")
                        .help("Taxon to be supressed, for modelling/training purposes")
                        .required(false)
                        .short("S")
                        .takes_value(true)
                        .long("supress_taxon"),
                        ),
        )
        .subcommand(
            SubCommand::with_name("batch_id")
                .about("classifies batch of samples reads")
                .version("0.2")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("bigsi")
                        .short("b")
                        .long("bigsi")
                        .required(true)
                        .takes_value(true)
                        .help("index to be used for search"),
                        )
                .arg(
                    Arg::with_name("query")
                        .help("tab-delimiteed file with samples to be classified [sample_name reads1 reads2 (optional)]")
                        .required(true)
                        .min_values(1)
                        .short("q")
                        .takes_value(true)
                        .long("query"),
                )
                .arg(
                    Arg::with_name("tag")
                        .help("tag to be include in output file names ")
                        .required(true)
                        .short("T") //
                        .takes_value(true)
                        .long("tag"),
                )
                .arg(
                    Arg::with_name("batch")
                        .help("Sets size of batch of reads to be processed in parallel (default 50,000)")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("batch"),
                )
                .arg(
                    Arg::with_name("threads")
                        .help("number of threads to use, if not set the maximum available number threads will be used")
                        .required(false)
                        .short("t")
                        .takes_value(true)
                        .long("threads"),
                )
                .arg(
                    Arg::with_name("down_sample")
                        .help("down-sample k-mers used for read classification, default 1; increases speed at cost of decreased sensitivity ")
                        .required(false)
                        .short("d")
                        .takes_value(true)
                        .long("down_sample"),
                )
                .arg(
                    Arg::with_name("fp_correct")
                        .help("Parameter to correct for false positives, default 3 (= 0.001), maybe increased for larger searches. Adjust for larger datasets")
                        .required(false)
                        .short("p")
                        .takes_value(true)
                        .long("fp_correct"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("kmers with nucleotides below this minimum phred score will be excluded from the analyses (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                )
                .arg(
                    Arg::with_name("bitvector_sample")
                        .help("Collects matches for subset of kmers indicated (default=3), using this subset to more rapidly find hits for the remainder of the kmers")
                        .required(false)
                        .short("B")
                        .takes_value(true)
                        .long("bitvector_sample"),
                        )
                .arg(
                    Arg::with_name("supress_taxon")
                        .help("Taxon to be supressed, for modelling/training purposes")
                        .required(false)
                        .short("S")
                        .takes_value(true)
                        .long("supress_taxon"),
                        ),
        )
        .subcommand(
            SubCommand::with_name("read_filter")
                .about("filters reads")
                .version("0.1")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("classification")
                        .short("c")
                        .long("classification")
                        .required(true)
                        .takes_value(true)
                        .help("tab delimited read classification file generated with the read_id subcommand"),
                        )
                .arg(
                    Arg::with_name("files")
                        .help("query file(-s)fastq.gz")
                        .required(true)
                        .min_values(1)
                        .short("f")
                        .takes_value(true)
                        .long("files"),
                )
                .arg(
                    Arg::with_name("taxon")
                        .help("taxon to be in- or excluded from the read file(-s)")
                        .required(true)
                        .short("t")
                        .takes_value(true)
                        .long("taxon"),
                        )
                .arg(
                    Arg::with_name("prefix")
                        .help("prefix for output file(-s)")
                        .required(true)
                        .short("p")
                        .takes_value(true)
                        .long("prefix"),
                        )
                .arg(
                    Arg::with_name("exclude")
                        .help("If set('-e or --exclude'), reads for which the classification contains the taxon name will be excluded")
                        .required(false)
                        .short("e")
                        .takes_value(false)
                        .long("exclude"),
                ),
          )
          .subcommand(
            SubCommand::with_name("merge")
                .about("merges (concatenates) indices")
                .version("0.1")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("index_1")
                        .short("1")
                        .long("index_1")
                        .required(true)
                        .takes_value(true)
                        .help("index to which index 2 will be concatenated"),
                        )
                .arg(
                    Arg::with_name("index_2")
                        .help("index to be concatenated to index 1")
                        .required(true)
                        .short("2")
                        .takes_value(true)
                        .long("index_2"),
                )
                .arg(
                    Arg::with_name("out_bigsi")
                        .short("b")
                        .long("out_bigsi")
                        .required(true)
                        .takes_value(true)
                        .help("name output index"),
                ),
          )
              .get_matches();
    if let Some(matches) = matches.subcommand_matches("build") {
        println!(" Ref_file : {}", matches.value_of("ref_file").unwrap());
        println!(" Bigsi file : {}", matches.value_of("bigsi").unwrap());
        println!("K-mer size: {}", matches.value_of("k-mer_size").unwrap());
        println!(
            "Bloom filter parameters: num hashes {}, filter size {}",
            matches.value_of("num_hashes").unwrap(),
            matches.value_of("length_bloom").unwrap()
        );
        let kmer = value_t!(matches, "k-mer_size", usize).unwrap_or(31);
        let bloom = value_t!(matches, "length_bloom", u64).unwrap_or(50_000_000);
        let hashes = value_t!(matches, "num_hashes", u64).unwrap_or(4);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let batch = value_t!(matches, "batch", usize).unwrap_or(300);
        let minimizer = value_t!(matches, "minimizer", usize).unwrap_or(0);
        let map = build::tab_to_map(matches.value_of("ref_file").unwrap().to_string());
        let filter = value_t!(matches, "filter", isize).unwrap_or(-1);
        let index = matches.value_of("bigsi").unwrap();
        if minimizer > 0 {
            println!("Build with minimizers, minimizer size: {}", minimizer);
            let (bigsi_map, colors_accession, n_ref_kmers) = build::build_multi_mini_bigvec(
                &map,
                bloom as usize,
                hashes,
                kmer,
                minimizer,
                threads,
                quality,
                filter,
                batch,
            );
            eprintln!("Saving index to file.");
            //fs::create_dir_all("/".to_string() + index).expect("could not initiate db");
            fs::DirBuilder::new()
                .recursive(true)
                .create("./".to_string() + index)
                .expect("could not initiate db");

            bigsi_map.save(&("./".to_string() + index + "/bigsi"));

            let parameters = build::Parameters {
                k_size: kmer,
                m_size: minimizer,
            };

            let serialized: Vec<u8> = serialize(&parameters).unwrap();
            let mut writer = fs::File::create(&("./".to_string() + index + "/parameters")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing serialized parameters for writing");

            let serialized: Vec<u8> = serialize(&colors_accession).unwrap();
            let mut writer = File::create(&("./".to_string() + index + "/colors")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing serialized colors for writing");

            let serialized: Vec<u8> = serialize(&n_ref_kmers).unwrap();
            let mut writer = File::create(&("./".to_string() + index + "/ref_kmers")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing kmer information for writing");
        //parameters
        /*};
        println!("Saving BIGSI to file.");
        let index = bigsi::BigsyMapMiniNew {
            map: bigsi_map,
            colors: colors_accession,
            n_ref_kmers: n_ref_kmers,
            bloom_size: bloom,
            num_hash: hashes,
            k_size: kmer,
            m_size: minimizer_value,
        };

        bigsi::save_bigsi_mini(
            &(matches.value_of("bigsi").unwrap().to_owned() + ".mxi"),
            &index,
        );*/
        /*
        bigsi::save_bigsi_mini(
            bigsi_map.to_owned(),
            colors_accession.to_owned(),
            n_ref_kmers.to_owned(),
            bloom,
            hashes,
            kmer,
            minimizer_value,
            &(matches.value_of("bigsi").unwrap().to_owned() + ".mxi"),
        );*/
        } else {
            let (bigsi_map, colors_accession, n_ref_kmers) = build::build_multi_bigvec_mutex(
                &map,
                bloom as usize,
                hashes,
                kmer,
                threads,
                quality,
                filter,
                batch,
            );
            eprintln!("Saving index to file.");
            //fs::create_dir_all("/".to_string() + index).expect("could not initiate db");
            fs::DirBuilder::new()
                .recursive(true)
                .create("./".to_string() + index)
                .expect("could not initiate db");

            bigsi_map.save(&("./".to_string() + index + "/bigsi"));

            let parameters = build::Parameters {
                k_size: kmer,
                m_size: minimizer,
            };

            let serialized: Vec<u8> = serialize(&parameters).unwrap();
            let mut writer = fs::File::create(&("./".to_string() + index + "/parameters")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing serialized parameters for writing");

            let serialized: Vec<u8> = serialize(&colors_accession).unwrap();
            let mut writer = File::create(&("./".to_string() + index + "/colors")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing serialized colors for writing");

            let serialized: Vec<u8> = serialize(&n_ref_kmers).unwrap();
            let mut writer = File::create(&("./".to_string() + index + "/ref_kmers")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing kmer information for writing");
        }
    }
    if let Some(matches) = matches.subcommand_matches("search") {
        let files1: Vec<_> = matches.values_of("query").unwrap().collect();
        let files2 = if matches.value_of("reverse").unwrap() == "none" {
            vec![]
        } else {
            matches.values_of("reverse").unwrap().collect()
        };
        //let files2: Vec<_> = matches.values_of("reverse").unwrap().collect();
        let filter = value_t!(matches, "filter", isize).unwrap_or(-1);
        let cov = value_t!(matches, "shared_kmers", f64).unwrap_or(0.35);
        let gene_search = matches.is_present("gene_search");
        let perfect_search = matches.is_present("perfect_search");
        let multi_fasta = matches.is_present("multi_fasta");
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let index = matches.value_of("bigsi").unwrap();
        //parameters
        let mut reader = BufReader::new(
            File::open(&(index.to_owned() + "/parameters")).expect("Can't open parameters!"),
        );
        let parameters: build::Parameters =
            deserialize_from(&mut reader).expect("can't deserialize");
        if parameters.m_size > 0 {
            eprintln!(
                "Error: An index with minimizers (.mxi) is used, but not available for this function"
            );
        } else {
            let bigsi_time = SystemTime::now();
            eprintln!("Loading index");
            let mut db = bigsi_rs::Bigsi::default();
            db.read(&(index.to_owned() + "/bigsi"));

            //parameters
            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/parameters")).expect("Can't open parameters!"),
            );
            let parameters: build::Parameters =
                deserialize_from(&mut reader).expect("can't deserialize");
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

            //let index = bigsi::read_bigsi(matches.value_of("bigsi").unwrap());
            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if perfect_search {
                if multi_fasta {
                    //make 'perfect' batch....
                    colorid::perfect_search::batch_search_mf(
                        files1,
                        &db,
                        &colors,
                        parameters.k_size,
                        cov,
                    )
                } else {
                    colorid::perfect_search::batch_search(
                        files1,
                        &db,
                        &colors,
                        parameters.k_size,
                        cov,
                    )
                }
            } else {
                if multi_fasta {
                    //make 'perfect' batch....
                    colorid::batch_search_pe::batch_search_mf(
                        files1,
                        &db,
                        &colors,
                        parameters.k_size,
                        cov,
                    )
                } else {
                    colorid::batch_search_pe::batch_search(
                        files1,
                        files2,
                        &db,
                        &colors,
                        &ref_kmers,
                        parameters.k_size,
                        filter,
                        cov,
                        gene_search,
                        quality,
                    )
                }
            }
        }
    }

    if let Some(matches) = matches.subcommand_matches("info") {
        let bigsi_time = SystemTime::now();
        eprintln!("Loading index");
        let index = matches.value_of("bigsi").unwrap();
        let mut db = bigsi_rs::Bigsi::default();
        db.read(&(index.to_owned() + "/bigsi"));
        //parameters
        let mut reader = BufReader::new(
            File::open(&(index.to_owned() + "/parameters")).expect("Can't open parameters!"),
        );
        let parameters: build::Parameters =
            deserialize_from(&mut reader).expect("can't deserialize");
        //colors
        let mut reader = BufReader::new(
            File::open(&(index.to_owned() + "/colors")).expect("Can't open parameters!"),
        );
        let colors: fnv::FnvHashMap<usize, String> =
            deserialize_from(&mut reader).expect("can't deserialize");
        //ref_kmers
        let mut reader = BufReader::new(
            File::open(&(index.to_owned() + "/ref_kmers")).expect("Can't open parameters!"),
        );
        let ref_kmers: fnv::FnvHashMap<String, usize> = deserialize_from(&mut reader).expect("can't deserialize");
        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index loaded in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }
        println!(
            "kmer size: {}\n minimizer size: {}",
            parameters.k_size, parameters.m_size
        );
        println!(
            "bloom size: {}\n number of hashes: {}",
            db.bigsi.len(),
            db.num_hashes
        );
        println!("accessions:");
        for (_key, value) in colors {
            let k_size = ref_kmers.get(&value).unwrap();
            println!(
                "{} {} {:.3}",
                value,
                k_size,
                false_prob(db.bigsi.len() as f64, db.num_hashes as f64, *k_size as f64)
            );
        }

        /*
        if index.ends_with(".mxi") {
            let bigsi = bigsi::read_bigsi_mini(index);
            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            println!(
            "BIGSI parameters:\nBloomfilter-size: {}\nNumber of hashes: {}\nK-mer size: {}\n minimizer size: {}\n",
            bigsi.bloom_size, bigsi.num_hash, bigsi.k_size, bigsi.m_size
        );
            println!("Number of accessions in index: {}", bigsi.colors.len());
            let mut accessions = Vec::new();
            for (_k, v) in bigsi.colors {
                accessions.push(v);
            }
            accessions.sort_by(|a, b| a.cmp(b));
            for a in accessions {
                let k_size = bigsi.n_ref_kmers.get(&a).unwrap();
                println!(
                    "{} {} {:.3}",
                    a,
                    k_size,
                    false_prob(
                        bigsi.bloom_size as f64,
                        bigsi.num_hash as f64,
                        *k_size as f64
                    )
                );
            }
        } else {
            let bigsi = bigsi::read_bigsi(index);
            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            println!(
                "BIGSI parameters:\nBloomfilter-size: {}\nNumber of hashes: {}\nK-mer size: {}",
                bigsi.bloom_size, bigsi.num_hash, bigsi.k_size
            );
            println!("Number of accessions in index: {}", bigsi.colors.len());
            let mut accessions = Vec::new();
            for (_k, v) in bigsi.colors {
                accessions.push(v);
            }
            accessions.sort_by(|a, b| a.cmp(b));
            for a in accessions {
                let k_size = bigsi.n_ref_kmers.get(&a).unwrap();
                println!(
                    "{} {} {:.3}",
                    a,
                    k_size,
                    false_prob(
                        bigsi.bloom_size as f64,
                        bigsi.num_hash as f64,
                        *k_size as f64
                    )
                );
            }
        }*/
    }
    if let Some(matches) = matches.subcommand_matches("read_id") {
        let bigsi_time = SystemTime::now();
        //let fq = matches.value_of("query").unwrap();
        let fq: Vec<_> = matches.values_of("query").unwrap().collect();
        let threads = value_t!(matches, "threads", usize).unwrap_or(0);
        let down_sample = value_t!(matches, "down_sample", usize).unwrap_or(1);
        let correct = value_t!(matches, "fp_correct", f64).unwrap_or(3.0);
        let fp_correct = 10f64.powf(-correct);
        let index = matches.value_of("bigsi").unwrap();
        let prefix = matches.value_of("prefix").unwrap();
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let batch = value_t!(matches, "batch", usize).unwrap_or(50000);
        let bitvector_sample = value_t!(matches, "bitvector_sample", usize).unwrap_or(3);
        let supress_taxon = matches.value_of("supress_taxon").unwrap_or("None");
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Can't initialize ThreadPoolBuilder");

        let mut db = bigsi_rs::Bigsi::default();
        db.read(&(index.to_owned() + "/bigsi"));

        //parameters
        let mut reader = BufReader::new(
            File::open(&(index.to_owned() + "/parameters")).expect("Can't open parameters!"),
        );
        let parameters: build::Parameters =
            deserialize_from(&mut reader).expect("can't deserialize");
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

        match bigsi_time.elapsed() {
            Ok(elapsed) => {
                eprintln!("Index loaded in {} seconds", elapsed.as_secs());
            }
            Err(e) => {
                // an error occurred!
                eprintln!("Error: {:?}", e);
            }
        }

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
                    prefix,
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
                    prefix,
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
                prefix,
                bitvector_sample,
                supress_taxon,
            );
        }
        colorid::reports::read_counts_five_fields(prefix.to_owned() + "_reads.txt", prefix);
    }
    if let Some(matches) = matches.subcommand_matches("batch_id") {
        let batch_samples = matches.value_of("query").unwrap();
        //let fq: Vec<_> = matches.values_of("query").unwrap().collect();
        let threads = value_t!(matches, "threads", usize).unwrap_or(0);
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Can't initialize ThreadPoolBuilder");
        let down_sample = value_t!(matches, "down_sample", usize).unwrap_or(1);
        let correct = value_t!(matches, "fp_correct", f64).unwrap_or(3.0);
        let fp_correct = 10f64.powf(-correct);
        let index = matches.value_of("bigsi").unwrap();
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let batch = value_t!(matches, "batch", usize).unwrap_or(50000);
        let bitvector_sample = value_t!(matches, "bitvector_sample", usize).unwrap_or(3);
        let tag = matches.value_of("tag").unwrap();
        let supress_taxon = matches.value_of("supress_taxon").unwrap_or("None");

        colorid::read_id_batch::read_id_batch(
            batch_samples,
            index,
            down_sample,
            fp_correct,
            batch,
            quality,
            bitvector_sample,
            tag,
            supress_taxon,
        );
    }
    if let Some(matches) = matches.subcommand_matches("read_filter") {
        let classification = matches.value_of("classification").unwrap();
        let files: Vec<_> = matches.values_of("files").unwrap().collect();
        let taxon = matches.value_of("taxon").unwrap();
        let prefix = matches.value_of("prefix").unwrap();
        let exclude = matches.is_present("exclude");
        let map = colorid::read_filter::tab_to_map(classification.to_string(), taxon);
        if files.len() == 1 {
            colorid::read_filter::read_filter_se(map, files, taxon, prefix, exclude);
        } else {
            colorid::read_filter::read_filter_pe(map, files, taxon, prefix, exclude);
        }
    }
    if let Some(matches) = matches.subcommand_matches("merge") {
        let index_1 = matches.value_of("index_1").unwrap();
        let index_2 = matches.value_of("index_2").unwrap();
        let index_out = matches.value_of("out_bigsi").unwrap();
        //1. read indices and parameters; check if parameters (# hashes, kmer, mmer size) are compatible
        //parameters
        let mut reader = BufReader::new(
            File::open(&(index_1.to_owned() + "/parameters")).expect("Can't open parameters!"),
        );
        let parameters_1: build::Parameters =
            deserialize_from(&mut reader).expect("can't deserialize");
        //colors
        let mut reader = BufReader::new(
            File::open(&(index_1.to_owned() + "/colors")).expect("Can't open parameters!"),
        );
        let mut colors_1: fnv::FnvHashMap<usize, String> =
            deserialize_from(&mut reader).expect("can't deserialize");
        //ref_kmers
        let mut reader = BufReader::new(
            File::open(&(index_1.to_owned() + "/ref_kmers")).expect("Can't open parameters!"),
        );
        let mut ref_kmers_1: fnv::FnvHashMap<String, usize> =
            deserialize_from(&mut reader).expect("can't deserialize");

        //parameters
        let mut reader = BufReader::new(
            File::open(&(index_2.to_owned() + "/parameters")).expect("Can't open parameters!"),
        );
        let parameters_2: build::Parameters =
            deserialize_from(&mut reader).expect("can't deserialize");

        //colors
        let mut reader = BufReader::new(
            File::open(&(index_2.to_owned() + "/colors")).expect("Can't open parameters!"),
        );
        let colors_2: fnv::FnvHashMap<usize, String> =
            deserialize_from(&mut reader).expect("can't deserialize");
        //ref_kmers
        let mut reader = BufReader::new(
            File::open(&(index_2.to_owned() + "/ref_kmers")).expect("Can't open parameters!"),
        );
        let ref_kmers_2: fnv::FnvHashMap<String, usize> =
            deserialize_from(&mut reader).expect("can't deserialize");
        //index_1
        let mut db_1 = bigsi_rs::Bigsi::default();
        db_1.read(&(index_1.to_owned() + "/bigsi"));
        //index_2
        let mut db_2 = bigsi_rs::Bigsi::default();
        db_2.read(&(index_2.to_owned() + "/bigsi"));
        if db_1.bigsi.len() != db_2.bigsi.len() {
            panic!("length indices not compatible (should be same)!");
        }
        if parameters_1 != parameters_2 {
            panic!("kmer/minimizer sizes not compatible!");
        }

        //2. concatenate bigsis
        db_1.merge(&db_2);
        println!("bit vector length after merge {}", db_1.bigsi[0].len());
        db_1.slim();
        //3. merge colors
        let accessions_1 = colors_1.len();
        for (key, value) in colors_2 {
            colors_1.insert(key + accessions_1, value);
        }
        for (key, value) in &colors_1 {
            eprintln!("{} {}", key, value);
        }
        //4. merge refkmers
        for (key, value) in ref_kmers_2 {
            ref_kmers_1.insert(key, value);
        }
        //5. safe new bigsi
        eprintln!("Saving index to file.");
        //fs::create_dir_all("/".to_string() + index).expect("could not initiate db");
        fs::DirBuilder::new()
            .recursive(true)
            .create("./".to_string() + index_out)
            .expect("could not initiate db");

        db_1.save(&("./".to_string() + index_out + "/bigsi"));

        let serialized: Vec<u8> = serialize(&parameters_1).unwrap();
        let mut writer = fs::File::create(&("./".to_string() + index_out + "/parameters")).unwrap();
        writer
            .write_all(&serialized)
            .expect("problems preparing serialized parameters for writing");

        let serialized: Vec<u8> = serialize(&colors_1).unwrap();
        let mut writer = File::create(&("./".to_string() + index_out + "/colors")).unwrap();
        writer
            .write_all(&serialized)
            .expect("problems preparing serialized colors for writing");

        let serialized: Vec<u8> = serialize(&ref_kmers_1).unwrap();
        let mut writer = File::create(&("./".to_string() + index_out + "/ref_kmers")).unwrap();
        writer
            .write_all(&serialized)
            .expect("problems preparing kmer information for writing");
    }
    Ok(())
}
