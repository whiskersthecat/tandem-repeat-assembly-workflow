# Peter Reifenstein, Alex Kozik
# Reference: https://github.com/whiskersthecat/tiger-paw

if [ "$#" -ne 2 ]; then
    echo "Usage: bash x-cut-transpose-x.sh variant_table.tab output_name"
    exit 2
fi

# -------------------------------

echo "cutting big table"
cut -f 1-7 --complement $1 > $1.Variants.$2.Columns

# -------------------------------

echo "transposing table"
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]"\t"
        for(i=2; i<=NR; i++){
            str=str""a[i,j];
        }
        print str
    }
}' $1.Variants.$2.Columns > $1.Variants.$2.Rows

# -------------------------------

echo "sorting"
sort -k2 $1.Variants.$2.Rows > $1.Variants.$2.Sorted

# -------------------------------

echo "grouping and counting"
cut -f 2 $1.Variants.$2.Sorted | sort | uniq -c | sort -k1,1nr > $1.Variants.$2.Uniq.Count