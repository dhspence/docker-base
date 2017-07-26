#!/usr/bin/perl

use strict;
use File::Basename;
use String::ShellQuote qw(shell_quote);

my $QCVERSION = basename($0,".pl");

my $BEDTOOLS="/usr/local/bin/bedtools";
my $SAMTOOLS="/usr/local/bin/samtools";
my $RSCRIPT="/usr/bin/Rscript";
my $COVERAGEPLOTSR="/usr/local/bin/CoveragePlots.R";

my $REFFASTA=$ARGV[0]; 
my $AMPLICONBED=$ARGV[1];
my $TARGETBED=$ARGV[2];
my $COVERAGEBED=$ARGV[3];
my $MINCOV1=$ARGV[4];
my $MINCOV2=$ARGV[5];
my $DEMUXFILE = $ARGV[6];
my $ABAM=$ARGV[7];
my $CBAM=$ARGV[8];
my $NAME=$ARGV[9];

open(F,">$NAME.qc.txt") || die;
print F "#sample: $NAME\t\n#qc version: $QCVERSION lowcov value: $MINCOV1 nocov value: $MINCOV2\n";
print F "#coverage bed file: $COVERAGEBED\n#amplicon bed file: $AMPLICONBED\n#target bed file: $TARGETBED\n";

my $flowcellreads = 0;
my $unknownreads = 0;
my $nsamples = 0;
open(D,$DEMUXFILE) || die;
while(<D>){
  chomp;
  my @l = split("\t",$_);
  $flowcellreads += $l[4];
  if(/unknown/){
    $unknownreads = $l[4];
  } else {
    $nsamples++;
  }
}

print F "FLOWCELL_READS\t$flowcellreads\n";
print F "FLOWCELL_SAMPLES\t$nsamples\n";
print F "UNKNOWN_READS\t$unknownreads\n";

my $TOTALREADS=`$SAMTOOLS view -c -F 0x900 $ABAM`; # no secondary or supp alignments
chomp $TOTALREADS;
my $ALIGNEDREADS=`$SAMTOOLS view -c -F 0x904 $ABAM`; # no secondary or supp alignments or unmapped reads
chomp $ALIGNEDREADS;
my $ONTARGETREADS=`$SAMTOOLS view -c -F 0x904 -L $TARGETBED $ABAM`; # overlap bed file
chomp $ONTARGETREADS;
my $CONDENSEDONTARGETREADS=`$SAMTOOLS view -c -F 0x904 -L $TARGETBED $CBAM`; # overlap bed file in condensed bam
chomp $CONDENSEDONTARGETREADS;

print F "TOTAL_READS\t$TOTALREADS\n";
print F "ALIGNED_READS\t$ALIGNEDREADS\n";
print F "PERCENT_ALIGNED_READS\t" . sprintf("%.1f",$ALIGNEDREADS / $TOTALREADS * 100) . "\n";
print F "ONTARGET_READS\t$ONTARGETREADS\n";
print F "PERCENT_ONTARGET_READS\t" . sprintf("%.1f",$ONTARGETREADS / $ALIGNEDREADS * 100) . "\n";
print F "CONDENSED_ONTARGET_READS\t$CONDENSEDONTARGETREADS\n";
print F "PERCENT_CONDENSED_ONTARGET_READS\t" . sprintf("%.1f",$CONDENSEDONTARGETREADS / $ALIGNEDREADS * 100) . "\n";

# histogram of the BC tag. anything over 25 reads is grouped together
open(BC,"$SAMTOOLS view -F 0x904 -L $TARGETBED $CBAM |") || die;
my %h = ();
while(<BC>){
    $_ =~ /BC:Z:(\d+)/;
    my $i=$1;
    $i=25 if $i>25;
    $h{$1}++;
}
close BC;
print F "READS_PER_UMI\t",join(",", map { $h{$_} || 0 } 1..25),"\n";

# get coverage depth
`$SAMTOOLS depth -d 100000 -a -q 20 -Q 20 -b $COVERAGEBED $CBAM | /usr/bin/awk -v OFS="\t" '{ print \$1,\$2-1,\$2,\$3; }' | $BEDTOOLS intersect -a stdin -b $COVERAGEBED -wo | cut -f 1-4,9-11 > $NAME.coverage.txt`;

# positions that are poorly covered (>=20, <50)
open(C,"/usr/bin/awk -v mincov1=$MINCOV1 -v mincov2=$MINCOV2 -v OFS=\"\t\" '\$4<mincov1 && \$4>=mincov2' $NAME.coverage.txt | $BEDTOOLS merge -i stdin -c 5,5,4,4 -o count_distinct,distinct,count,collapse | /usr/bin/awk -v OFS=\"\t\" '{ print \"LOWCOV\",\$0; }' |") || die;
while(<C>){ 
    print F $_;
}
close C;

# positions that are basically not covered (<20)
open(C,"/usr/bin/awk -v mincov2=$MINCOV2 -v OFS=\"\t\" '\$4<mincov2' $NAME.coverage.txt | $BEDTOOLS merge -i stdin -c 5,5,4,4 -o count_distinct,distinct,count,collapse | /usr/bin/awk -v OFS=\"\t\" '{ print \"NOCOV\",\$0; }' |") || die;
while(<C>){ 
    print F $_;
}
close C;

# amplicon length and gc
my $cmd = "$BEDTOOLS nuc -fi $REFFASTA -bed $AMPLICONBED > nucs.bed && $SAMTOOLS view -uf 0x2 -F 0x900 $CBAM | $SAMTOOLS sort -n - | $BEDTOOLS bamtobed -mate1 -bedpe -i stdin 2> /dev/null | /usr/bin/awk -v OFS=\"\\t\" '{ if (\$9==\"+\"){ print \$1,\$2,\$6,\$7,\$8,\$9; } else { print \$1,\$5,\$3,\$7,\$8,\$9; } }' | $BEDTOOLS sort -i stdin | $BEDTOOLS intersect -r -f .50 -s -a nucs.bed -b stdin -wo | /usr/bin/awk -v OFS=\"\\t\" '{ print \$0,sqrt((\$17-\$2)^2)+sqrt((\$18-\$3)^2); }' | sort -k 19,19 -k 23n,23 | $BEDTOOLS groupby -g 19 -full -c 23 -o first | /usr/bin/awk '\$NF<20' | sort -k 1,1 -k 2n,2 -k 4,4 | $BEDTOOLS groupby -i stdin -g 1,2,3,4 -c 4,8 -o count,first > $NAME.amplicon_counts.txt && rm -f nucs.bed";
`$cmd`;

# make plots, return amplicon-level coverage by length and GC and gene-level coverage information
open(C,"$RSCRIPT $COVERAGEPLOTSR $NAME |") || die; #2> /dev/null # && rm -f amplicon_counts.txt coverage.txt
while(<C>){ 
    print F $_; 
}
close C;

close F;
