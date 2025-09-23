# convert the VCF (e.g. from freebayes) to a tab delimited format like from CLC

import sys

header = "\"Reference Position\"	\"Type\"	\"Length\"	\"Reference\"	\"Allele\"	\"Linkage\"	\"Zygosity\"	\"Count\"	\"Coverage\"	\"Frequency\"	\"Forward/reverse balance\"	\"Average quality\"	\"Overlapping annotations\"	\"Coding region change\"	\"Amino acid change\""

if(len(sys.argv) < 3):
    print("Usage: Python3 variants_to_vcf.py variants.var reference_name")
    exit()
    
reference_name = sys.argv[2]
variants_file = open(sys.argv[1],'r')
line = variants_file.readline()
vcf_file = open(sys.argv[1] + ".VCF",'w')

vcf_file.write("##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n")

n = 0
for line in variants_file.readlines():
    n += 1
    line = line.split()
    vcf_file.write(reference_name + "\t" + line[0] + "\tPOL" + str(n) + "\t" + line[3] + "\t" + line[4] + "\t.\tPASS\t.\n")
