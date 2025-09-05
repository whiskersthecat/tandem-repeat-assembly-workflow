# 00_requirements


To run *traw*, you will need two large DNA sequencing datasets, ideally both with coverage at least 50x. For *long.fa*, the more the better.

* ***short.fa***. Reads are up to 30,000 base pairs and each spans at least one segment (see definitions), and are very accurate. E.g. Pacbio HiFi Reads. They are assumed to contain no chimeras or duplexes.

* ***long.fa***. E.g. Reads are longer than 30,000 base pairs and can span at least three segments, ideally up to 20 or more. They are assumed to be innacurate, contain duplexes, and may contain chimeras. Oxford Nanopore Reads. 

These **definitions** are used throughout the workflow.

* ***segment***: one copy of the tandem repeat. A single long read should span many segments.
* ***chimera***: a single read in which two different (unrelated) molecules are concatenated together.
* ***duplex***: a single read in which the same molecule is read twice, with forward and reverse strand concatenated together.

### 00_A. Prepare reads

Note: in all commands presented in the workflow, text in \* \* should be replaced with user parameters.

```bash
mkdir reads
cp *path_to_short_reads* reads/short.fa
cp *path_to_long_reads* reads/long.fa
```


# 01_create-reference-segment

```bash
mkdir reference
```


# 02_prepare-reads

Goal: generate new collections of reads which contain 


### 02_A. Collect all reads with blast hits to your reference segment.

```bash
mkdir blastN
makeblastdb -in reads/short.fa -out reads/short.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
makeblastdb -in reads/long.fa -out reads/long.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
```
```bash
mkdir blastN
blastn -task blastn -query reference/reference.fa -db reads/short.BDB -out blastN/reference.vs.short.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads 100
blastn -task blastn -query reference/reference.fa -db reads/long.BDB -out blastN/reference.vs.long.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads 100
```

Collect large blast hits and reformat into table.

```bash
awk -F'\t' '$6 >= 2000' blastN/reference.vs.short.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.short.blastN.large
awk -F'\t' '$6 >= 2000' blastN/reference.vs.long.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.long.blastN.large
```

```bash
python3 02_chopDuplexAndOrientReads.py reads/short.fa blastN/reference.vs.short.blastN.large
python3 02_chopDuplexAndOrientReads.py reads/long.fa blastN/reference.vs.long.blastN.large
```




### 02_B. Orient all of reads in same direction, and check for duplex reads.

Note: Some tandem repeat clusters have orientation switches, and some do not. Distinguishing real orientation switches from duplex reads is non-trivial. 






# 03_read-segmentation

