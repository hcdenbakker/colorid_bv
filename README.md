# colorid_bv
Colorid using a new BIGSI implementation

# MLST

## Create the database from assemblies

Set some variables to make it easier downstream

```bash
fofn=$(realpath ./assemblies.fofn)
db=colorid.31.bxi
cgMLST=$HOME/PasteurcgMLST/*.fasta
```

Get the file of filenames and then build the database

```bash
\ls files/*.fasta | xargs -n 1 bash -c 'b=$(basename $0 .fasta); path=$(realpath $0); echo -e "$b\t$path";' > $fofn
colorid_bv build -b $db -s 30000000 -n 2 -k 31 -t 12 -r $fofn
```

## query

```bash
colorid_bv search -b $db -q $cgMLST -ms > alleles.txt
```

## post processing

Generates a bunch of reports and allele files, including the percentage of alleles not called or called with multiple alleles.

```bash
process_MLST.py alleles.txt colorid_
```
