#!/bin/bash

QCVERSION=$(basename $0 .sh)

BEDTOOLS="/usr/local/bin/bedtools"
SAMTOOLS="/usr/local/bin/samtools"
RSCRIPT="/usr/bin/Rscript"
COVERAGEPLOTSR=$(dirname $0)"/CoveragePlots.R"

REFFASTA=$1 
AMPLICONBED=$2
TARGETBED=$3
COVERAGEBED=$4
MINCOV1=$5
MINCOV2=$6
ABAM=$7
CBAM=$8
NAME=$9

echo -e "#sample: $NAME\t\n#qc version: $QCVERSION lowcov value: $MINCOV1 nocov value: $MINCOV2"
echo -e "#coverage bed file: $COVERAGEBED\n#amplicon bed file: $AMPLICONBED\n#target bed file: $TARGETBED"

TOTALREADS=$($SAMTOOLS view -c -F 0x900 $ABAM) # no secondary or supp alignments
ALIGNEDREADS=$($SAMTOOLS view -c -F 0x904 $ABAM) # no secondary or supp alignments or unmapped reads
ONTARGETREADS=$($SAMTOOLS view -c -F 0x904 -L $TARGETBED $ABAM) # overlap bed file
CONDENSEDONTARGETREADS=$($SAMTOOLS view -c -F 0x904 -L $TARGETBED $CBAM) # overlap bed file in condensed bam

echo -e "TOTAL_READS\t$TOTALREADS"
echo -e "ALIGNED_READS\t$ALIGNEDREADS"
echo -e "PERCENT_ALIGNED_READS\t"$(echo "$ALIGNEDREADS / $TOTALREADS * 100" | bc -l)
echo -e "ONTARGET_READS\t$ONTARGETREADS"
echo -e "PERCENT_ONTARGET_READS\t"$(echo "$ONTARGETREADS / $ALIGNEDREADS * 100" | bc -l)
echo -e "CONDENSED_ONTARGET_READS\t$CONDENSEDONTARGETREADS"
echo -e "PERCENT_CONDENSED_ONTARGET_READS\t"$(echo "$CONDENSEDONTARGETREADS / $ALIGNEDREADS * 100" | bc -l)

# histogram of the BC tag. anything over 25 reads is grouped together
$SAMTOOLS view -F 0x904 -L $TARGETBED $CBAM | /usr/bin/perl -ne '{/BC:Z:(\d+)/; $i=$1; $i=25 if $i>25; $h{$1}++;} END { print "READS_PER_UMI\t",join(",", map { $h{$_} || 0 } 1..25),"\n"; }'

# median coverage by gene and bw plot
$SAMTOOLS depth -d 100000 -a -q 20 -Q 20 -b $COVERAGEBED $CBAM | /usr/bin/awk -v OFS="\t" '{ print $1,$2-1,$2,1+c++,$3; }' | $BEDTOOLS intersect -a $COVERAGEBED -b stdin -wo | cut -f 1-7,11,12 > coverage.txt 

# positions that are poorly covered (>=20, <50)
/usr/bin/awk -v mincov1=$MINCOV1 -v mincov2=$MINCOV2 -v OFS="\t" '$9<mincov1 && $9>=mincov2 { $2=$2+$8-1; $3=$2+$8; print; }' coverage.txt | $BEDTOOLS merge -i stdin -c 5,5,9,9 -o count_distinct,distinct,count,collapse | /usr/bin/awk -v OFS="\t" '{ print "LOWCOV",$0; }'

# positions that are basically not covered (<20)
/usr/bin/awk -v mincov2=$MINCOV2 -v OFS="\t" '$9<mincov2 { $2=$2+$8-1; $3=$2+$8; print; }' coverage.txt | $BEDTOOLS merge -i stdin -c 5,5,9,9 -o count_distinct,distinct,count,collapse | /usr/bin/awk -v OFS="\t" '{ print "NOCOV",$0; }'

# amplicon length and gc
$BEDTOOLS intersect -r -f .50 -s -a <($BEDTOOLS nuc -fi $REFFASTA -bed $AMPLICONBED) -b <($SAMTOOLS view -uf 0x2 -F 0x900 $CBAM | $SAMTOOLS sort -n - | $BEDTOOLS bamtobed -mate1 -bedpe -i stdin 2> /dev/null | /usr/bin/awk -v OFS="\t" '{ if ($9=="+"){ print $1,$2,$6,$7,$8,$9; } else { print $1,$5,$3,$7,$8,$9; } }' | $BEDTOOLS sort -i stdin) -wo | /usr/bin/awk -v OFS="\t" '{ print $0,sqrt(($17-$2)^2)+sqrt(($18-$3)^2); }' | sort -k 19,19 -k 23n,23 | $BEDTOOLS groupby -g 19 -full -c 23 -o first | /usr/bin/awk '$NF<20' | sort -k 1,1 -k 2n,2 -k 4,4 | $BEDTOOLS groupby -i stdin -g 1,2,3,4 -c 4,8 -o count,first > amplicon_counts.txt

# make plots, return amplicon-level coverage by length and GC and gene-level coverage information
$RSCRIPT $COVERAGEPLOTSR $NAME #2> /dev/null # && rm -f amplicon_counts.txt coverage.txt
