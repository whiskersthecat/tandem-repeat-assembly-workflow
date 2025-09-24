# 0. introduction


These **definitions** are used throughout the workflow.

* ***segment***: one repeat copy of a tandem repeat. A single long read should span many segments.
* ***chimera***: a single read in which two different (unrelated) molecules are concatenated together.
* ***duplex***: a single read in which the same molecule is read twice, with forward and reverse strand concatenated together.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>




To run **tiger-paw**, you will need two large DNA sequencing datasets, ideally both with coverage at least 50x. For *long.fa*, the more the better.

* ***short.fa***. Reads are up to 30,000 base pairs and each spans at least one segment, and are **very accurate**. They are assumed to contain no chimeras or duplexes. For example, *Pacbio HiFi Reads*. 
> These reads will be used to generate common Haplotypes that will be used in assembly.


* ***long.fa***. Reads are longer than 30,000 base pairs and can span at **least three segments, ideally up to 20 or more**. They are assumed to be innacurate, contain duplexes, and may contain chimeras. For example, *Oxford Nanopore Reads*. 
> These reads will be used to actually assemble the tandem repeat cluster. 

Note: in the case that you have multiple datasets of either type, it is best to keep them seperate and run them in parallel. This helps confirm phenomena that is replicated among all datasets.


<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Dependencies:
- minimap2 (https://github.com/lh3/minimap2)
- blastn Version 2.16.0 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)
   - Note: different blastn Versions have very different performance. Newest blast version 2.17.0 struggled in identifying duplexes in HiFi reads.
- Python3
- linux

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 0A. Copy reads

Note: in all commands presented in the workflow, text in \* \* should be replaced with user parameters.

```bash
mkdir reads
cp *path_to_short_reads* reads/short.fa
cp *path_to_long_reads* reads/long.fa
```


# 1. creating reference segment

```bash
mkdir reference
```


# 2. preparing reads

Overview: for each dataset, generate a refined subset of reads which contain repeat segments.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 2a. Find reads with BLAST hits to reference segment.

```bash
mkdir blastN
mkdir blastN/DB
```
```bash
makeblastdb -in reads/short.fa -out blastN/DB/short.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
makeblastdb -in reads/long.fa -out blastN/DB/long.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
```
Modify \*THREADS\* based on your available computing resources.
```bash
blastn -task blastn -query reference/reference.fa -db blastN/DB/short.BDB -out blastN/reference.vs.short.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
blastn -task blastn -query reference/reference.fa -db blastN/DB/long.BDB -out blastN/reference.vs.long.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 2b. Reformat BLAST hits.

Collect large blast hits and reformat table.

For ***short*** dataset, the value \*SIZE\* should be adjusted to be roughly 60% of your segment length.
For ***long*** dataset, the value \*SIZE\* should be adjusted to be roughly 20% of your segment length.\

```bash
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.short.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.short.sorted.blastN
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.long.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.long.sorted.blastN
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 2c. Orienting of reads and finding duplexes.

In tandem repeats, segments may switch orientation at certain loci.
Notes:

- If this phenomenon is not observed in ***short*** dataset, this phenomenon does not exist. 
In this case, any switches in orientation observed in the ***long*** dataset represent duplex reads.
These duplex reads should be chopped, and the full length read should be removed from the dataset.

- On the other hand, if such a phenomenon is observed in the ***short*** dataset, these likely truly exist. There are two cases:
   + The orientation switch occurs in a single or few segments. In this case, duplex reads from the ***long*** dataset should be chopped like above.
   + The orientation switch occurs systematically at a certain breakpoint. In this case, distinguishing duplex reads in the long.fa dataset from real reads spanning this region is not possible at this point. Reads marked as duplex in ***long*** dataset should be used with caution during assembly.

Run the following script on both datasets.
```bash
python3 02_chopDuplexAndOrientReads.py reads/short.fa blastN/reference.vs.short.sorted.blastN --include_original False
python3 02_chopDuplexAndOrientReads.py reads/long.fa blastN/reference.vs.long.sorted.blastN --include_original True
```

This script will generate multiple output files.
1. reads/short.fa.**report** a count of how many reads were in each of 5 categories
2. reads/short.fa.**readids** the categorization of each read
3. reads/short.fa.**chop.oriented.fa** an output fasta file that includes modified reads:
   - With --include_original True (the default), any reads marked as duplexes will have two corresponding reads in the output file:
      - read_name: the modified read.
      - read_name_dp (duplex perfect) or read_name_dbf (duplex bias forward) or read_name_dbr (duplex bias reverse). This read has the sequence of the original read, and indicates its duplex categorization.
   - With --include_original False, only modified reads are output

# 3. segmenting reads

Overview: for each dataset, cut reads containing tandem repeats into segments.\
Consider this simplified example with a single 5bp reference segment of 'AAATG':

**input_reads.fa:**
```tsv
read_01  G AAATG AAATG AATTG AAATG AAATG AA
read_02  CCCCGCCC AAATG A
```
**blastHits.blastN:**
```tsv
read_01  1 2---- 3---- 4---- 5---- 6---- 7-
read_02           1----
```
**output_reads.fa:**
```tsv
read_01_segment_02 AAATG
read_01_segment_03 AAATG
read_01_segment_04 AATTG
read_01_segment_05 AAATG
read_01_segment_06 AAATG
read_01_segment_07 AA

read_02_segment_01 AAATG
```

Notes:
- the first segment of read_01 was removed because it was small (only 1bp)
- read_02 might contain the linkage between the chromosome and the tandem repeat. Alternatively, it could be a chimeric read.
- read_01_segment_04 contains a variant (T instead of A). This read could be a variant and serve useful in assembly. Alternatively, it could be a sequencing error that will confound downstream assembly.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 3a. Redo blast hits on your new dataset with oriented and chopped reads.


```bash
makeblastdb -in reads/short.fa.chop.oriented.fa -out blastN/DB/short.fa.chop.oriented.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
makeblastdb -in reads/long.fa.chop.oriented.fa -out blastN/DB/long.chop.oriented.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
```
Modify \*THREADS\* based on your available computing resources.
```bash
blastn -task blastn -query reference/reference.fa -db blastN/DB/short.fa.chop.oriented.BDB -out blastN/reference.vs.short.fa.chop.oriented.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
blastn -task blastn -query reference/reference.fa -db blastN/DB/long.chop.oriented.BDB -out blastN/reference.vs.long.fa.chop.oriented.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 3b. Reformat BLAST hits.

Like step 2b, collect large blast hits and reformat table.

For both datasets, dataset, the value \*SIZE\* should be adjusted to be roughly 20% of your reference segment length.\
```bash
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.short.fa.chop.oriented.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.short.fa.chop.oriented.sorted.blastN
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.long.fa.chop.oriented.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.long.fa.chop.oriented.sorted.blastN
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 3c. Segment the reads based on the BLAST hits.

Run the following script on both datasets.\

The value for \*SIZE\* determines the tolerance to include the segments at the edges of your reads.\
A larger (stricter) value will produce fewer segments.\
A smaller value will produce more partial (truncated) terminal segments.\

The recommended value for \*SIZE\* is different for your two datasets:
- For the short dataset, it should be very strict. At least 90% of your reference segment length (e.g. -l 9000 for a 10kb reference segment).
- For the long dataset, it can be more lenient. At to be least 70% of your reference segment length.

```bash
python3 03_extractSegments.py reads/short.fa.chopped.oriented.fa blastN/reference.vs.short.fa.chop.oriented.sorted.blastN -l *SIZE*
python3 03_extractSegments.py reads/long.fa.chopped.oriented.fa blastN/reference.vs.long.fa.chop.oriented.sorted.blastN -l *SIZE*
```

This script will generate multiple output files.
1. reads/short.fa.chopped.oriented.fa.**report** a quick summary
2. reads/short.fa.chopped.oriented.fa.**segmented.fa** an output fasta file that contains sequences of reads chopped into segments
3. reads/short.fa.chopped.oriented.fa.**undetected.fa** an output fasta file that includes sequences occuring between blast hits in a single read, that are of a minimum size
   - Small gaps between blast hits are common and not neccesarily indicative of a unique sequence (that is not a repeat segment)
   - The parameter --undetectedSegmentMinimumLength can be adjusted to count more or less of these as real features
4. reads/short.fa.chopped.oriented.fa.**allregions.fa** is a combination of file 2 and 3


# 4. generating haplotypes

Overview:
By identifying high frequency variant sites in the reference segment, types of segments with specific variants, called a ***haplotype***, can be identified.\

Using these, each segment in error-prone **long** reads can categorized into one of haplotypes.\

In essence, this washes away every other base pair in the segment, allowing focus on true biological variation between repeats and ignoring almost all sequencing error in every read.


<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4a. Prepare folders and rename segments file.
```bash
mkdir alignments
mv reads/short.fa.chop.oriented.fa.segmented.fa reads/short.segments.fa
mv reads/long.fa.chop.oriented.fa.segmented.fa reads/long.segments.fa
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4b. Align all of the segments both datasets to the reference segment.

```bash
minimap2 -t *THREADS* --eqx -a reference/reference.fa reads/short.segments.fa > alignments/short.segments.sam
minimap2 -t *THREADS* --eqx -a reference/reference.fa reads/long.segments.fa > alignments/long.segments.sam
```
Note: you won't need the alignment for the long.segments until part 5.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4c. Call variants on the alignment of the ***short.segments.sam*** alignment.

```bash
python3 04_highCoverageVariantCaller.py alignments/short.segments.sam reference/reference.fa --MODE SNP --min_frequency 0.02
```

Experiment with different min_frequency parameters for calling variants. \
It is recommended to set the min_frequency so that at least 10 biological segments in the tandem repeat cluster contain this variant.

You can estimate the total segments in your tandem repeat cluster.
```tsv
TOTAL_SEGMENTS_IN_TRC = total_segments_in_short_dataset / genomic_coverage_of_short_dataset
```
Then you can set min_frequency to be
```tsv
min_frequency = 10 / TOTAL_SEGMENTS_IN_TRC
```
For example, if TOTAL_SEGMENTS_IN_TRC = 1000 segments, this implies a --min_frequency of 0.01 requires at least 10 segments have this variant.



You may also call insertion and deletion variants. You may rename the output file to prevent overwriting.\
```bash
python3 04_highCoverageVariantCaller.py alignments/short.segments.sam reference/reference.fa --MODE INDEL --min_frequency 0.02
```

Notes:
- Though insertions and deletion variants are common in repeat clusters, they are difficult to precisely measure in noisy long reads.
- Therefore, it is recommended to use simple SNV and MNV variants, with the option --MODE SNP.


It is also highly recommended to visualize variants to ensure that they make sense.\
If there are variants in homopolymer regions, they are probably not reliable to use (even if they may be real). You may opt to duplicate the output file and remove variants deemed this way.


<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4d. Find nucleotides at important locations in alignments.

Consider that your variants file from 4c contains these four variants:
```tsv
Reference Position	Type	Length	Reference	Allele	Count	Coverage	Frequency
500	SNV	1	C	T	678	30014	2.2589458252
1000	SNV	1	T	G	1305	30014	4.3479709468
1500	SNV	1	C	T	1313	30014	4.3746251749
2000	MNV	2	TC	GG	1301	30014	4.3346438328
```
For any given segment, the ***trace*** of the segment is the nucleotides at these 4 variant positions.

For example, based on the frequencies above, this is a good prediction of the most common trace among the 30014 segments:
```tsv
CTCTC
```
Notes:
- Although all these variant occur at very low frequency, importantly they are well within the boundaries of not being systematic sequencing errors. 
Therefore, they represent real, rare, biological polymorphism in this tandem repeat cluster.
- Any MNV variants will take up multiple nucleotides in the trace, and are sorted to the end of the string.


Run this script.
```bash
python3 04_findVariants.py alignments/short.segments.sam alignments/short.segments.sam.VAR.tab reference/reference.fa
```
This will output a table which contains the nucleotide at each of the variant locations specified by the VAR file.

Next run this script.
```bash
bash 04_findTrace.sh alignments/short.segments.sam.variant_table.tab var1
```
Note:
- You can adjust **var1** to be another value for a different set of variants

This script will assemble the trace of each segment, and count the occurences of each trace.\

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4e. Creating haplotypes

```bash
mkdir haplotypes
```

Decide the minimum number of segments you wish to support one trace for it to be real (MIN_SUPPORT) and considered as a haplotype.\
It is recommended to choose \*MIN_SUPPORT\* such that:
```tsv
MIN_SUPPORT >= genomic_coverage_of_short_dataset
```

```bash
awk '$1 >= *MIN_SUPPORT* {printf "HAPLOTYPE_%03dC%04d\t%s\n", NR, $1, $2 }' alignments/short.segments.sam.VAR.tab.Traces.Counts > haplotypes/short.segments.sam.VAR.tab.Traces.Counts.Named
```

Its very useful to rename the haplotypes based on their phylogenetic clustering.\
First convert to fasta format.
```bash
awk '{print ">" $1 "\n" $2}' alignments/short.segments.sam.VAR.tab.Traces.Counts.Named > alignments/short.segments.sam.VAR.tab.Traces.Counts.Named.fa
```
Now run this program and input your fasta file:
```tsv
https://www.ebi.ac.uk/jdispatcher/phylogeny/simple_phylogeny
```
Download the output as ".tree" format into the alignments .\
Rename the file to "haplotypes.var1.tree"\
Run this script to map the original haplotype names to new haplotype names based on their clustering group.
```bash
python3 04_phylogenyLabeler.py alignments/haplotypes.var1.tree
```

```bash
mkdir haplotypes
sort alignments/haplotypes.var1.tree.RENAMED.tsv > alignments/haplotypes.var1.tree.RENAMED.sorted.tsv
join alignments/haplotypes.var1.tree.RENAMED.sorted.tsv alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces.Counts.Named | awk '{print $2"\t"$4}' > haplotypes/haplotypes.var1.tsv
```

The haplotypes file should look something like this:
```tsv
1A	CTCGCTAGTCAGAGTGAGTTCCTTAGATGACTTGGGGCCTCGAGATT
1D	CTCGCTAGTCAGAGTGAGTCCCTTAGATGACCTGTAGCCTCGAGATT
1B	CTCGCTAGTCAGAGTGAGTTCCTTAGATGACTTGGG---TCGAGATT
1A	CTCGCTAGTCAGAGTGAGTTCCTTAGATGACTTGGGG--TCGAGATT
7L	CTCGCTAGTCAGAGTGAGTCCCTTAGATGACCTGTA---TCGAGATT
6D	CGTAACGAGACAGACAGTCTGTACGGCTAGTCCATGTCCGGTCTCAA
7I	CTCGCTAGTCAGAGTGAGTCCCTTAGATGACCTGTAG--TCGAGATT
...
```

Notes:
- It is highly recommended to make two distinct set of haplotypes, at different depths of variant frequency.
- For example, one at 0.01 frequency, and one at 0.02 frequency.
- Higher frequencies will produce tighter haplotypes with smaller traces, that are easier to categorize and more robust to sequencing errors.
- Lower frequencies will allow for rarer variations to be used, at the cost of decreased accuracy in haplotype assigning (step 5).


# 5. creating blocks

A ***block*** is a single read listed as the haplotypes of its constituent segments, e.g.

```tsv
read_01  1A 1A 1D 1A 1A
read_03  1B 1D 
read_04  3A
```

Notes:
- 

### 5a. Assigning haplotypes to each segment

Generate the trace of each segment in the long.fa dataset.

```bash
python3 04_findVariants.py alignments/long.segments.sam alignments/long.segments.sam.VAR.tab reference/reference.fa
bash 04_findTrace.sh alignments/long.segments.sam.variant_table.tab var1
```

