#!/bin/bash

# collect primary contigs that have haloptigs assembled
contigs=`grep ">" $1 | cut -f2 -d ">" | cut -f1 -d"_" | sort | uniq -c | awk '$1>1 {print $2}'`

for contig in $contigs; do
    python3 extract_fasta.py ${contig} $1
    nucmer -prefix=$contig ${contig}_p.fa ${contig}_h.fa
    delta-filter -r -1 ${contig}.delta | show-coords -r -T -l -d -c /dev/stdin | awk 'NR>4' | sort -k14,14 -k15,15 -k1n,1 >  ${contig}.coords
    python3 aggressive_syntany.py ${contig}.coords ${contig} > ${contig}.loc.coords
done

cat *.loc.coords > loc.haplotig.coords
