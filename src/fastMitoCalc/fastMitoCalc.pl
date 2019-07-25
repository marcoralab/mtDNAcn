#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

my $chr = 0;
my $N = 3000;
my $region_size = 1000;
my $bedfile;
my $S;
my $prefix = "";
my $mt = "MT";

my @chr_sizes = (249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                 146364022,141213431,135534747,135006516,133851895,115169878,107349540,
                 102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566);
my $TotalChromLength = 2881033286;

sub usage() {
   print STDERR<<"   EOF";

   usage: perl $0 [options] -f bamfile -w workdir -p path_to_BaseCoverage
   Note: This version was written based on the GRCh37/hg19 chromosome sizes. 

     -f file     : bam file. Required
     -w work_dir      : working directory. Required
     -p BaseCoverage: BaseCoverage(complete path). Required
     -c chrID    : Estimate mt copy Number based on one chromosome chrID. No Random sampling. Optional.
     -n N        : Number of regions. default: 3000
     -s region_size  : sampled region size in random sampling. default: 1000.
     -b bedfile  : user supplied autosomal coordinates in bed format.
     -i          : keep the randomly generated coordinate file (bed) for checking. Optional
     -d          : keep the  depth file in the selected area. Optional
     -e prefix   : prefix of chromosome names (e.g. chr if names are like chr1,chr2,chr22) if chromesomes names are not 1,2,..22,MT
     -m mtname   : Mitochondria name if not MT (e.g, use -m chrM). 
     -g 37 or 38   : Human Genome version used in the bam (37 for GRCh37/hg19 or 38 for GRCh38/hg38). Default 37.
     -h          : this (help) message

   example: $0 -f /home/user/genome/bam_files/100.bam -w workDir -i -p /home/user/genome/bin/BaseCoverage

   EOF
   exit;
}


my %options=();
my $opt_string = 'hdic:f:w:n:s:b:p:e:m:g:';
getopts($opt_string, \%options ) or usage();

usage() if (defined $options{h});

if(!$options{f} || !$options{w} || !$options{p}) {
   print "Note: -f -w and -p are required parameters\n";
   usage();
}

my $bamfile = $options{f};
my $workDir = $options{w};
my $prog = $options{p};

my $bai = "$bamfile.bai";
my $bai2 = $bamfile;
substr($bai2,-1,1) = "i";

if(! -f $bai && ! -f $bai2) {
  print "error: You should have index file named as\n";
  print "\t$bai\n";
  print "or\n";
  print "\t$bai2\n\n";
  exit;
}

if($options{g}) {
  my $grch = $options{g};
  if($grch!=37 && $grch!=38) {
     print "use -g 37 for GRCh37 or -g 38 for GRCh38\n";
     exit;
  }
  if($grch==38) {
    @chr_sizes = (248956422,242193529,198295559,190214555,181538259,170805979,159345973,
                 145138636,138394717,133797422,135086622,133275309,114364328,107043718,
                 101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468);
    $TotalChromLength = 2875001522;
  }
}
if($options{s}) {
  $region_size = $options{s};
}
if($options{c}) {
  $chr = $options{c};
}
if($options{n}) {
  $N = $options{n};
}
if($options{s}) {
  $S = $options{s};
}
if($options{b}) {
  $bedfile = $options{b};
  $N =`wc -l < $bedfile`;
}
if($options{e}) {
  $prefix = $options{e};
}
if($options{m}) {
  $mt = $options{m};
}

my (@bam_path) = split(/\//,$bamfile);
my ($bamID) = split(/\./,$bam_path[$#bam_path]);

`mkdir -p $workDir`;

if($options{c}) {
  if( ($chr=~ /^\d+?$/) && $chr>=1 && $chr<=22) {
    `samtools view -b -q 20 -F 0x0704 $bamfile -o $workDir/$bamID.small.bam $prefix$chr $mt`;
  }
  elsif($chr=~ /\D/) {
    `samtools view -b -q 20 -F 0x0704 $bamfile -o $workDir/$bamID.small.bam $chr $mt`;
  }
  else {
     print "-c must be followed by the chromosome number, such as chr10, or 10\n";
     exit;
  }
}
elsif($options{b}) {
  convert_bed($bedfile);
  if($N<=30000) {
     my $str=make_str($bedfile);
     `samtools view -b -q 20 -F 0x0704 $bamfile $str -o $workDir/$bamID.small.bam`;
  }
  else {
    `samtools view -b -q 20 -F 0x0704 $bamfile -o $workDir/$bamID.small.bam -L $workDir/$bamID.bed`;
  }
}
else {
  my $str = make_bed();
  if($N<=30000) {
   `samtools view -b -q 20 -F 0x0704 $bamfile $str -o $workDir/$bamID.small.bam`;
  }
  else {
    `samtools view -b -q 20 -F 0x0704 $bamfile -o $workDir/$bamID.small.bam -L $workDir/$bamID.bed`;
  }
}

`samtools sort -O BAM $workDir/$bamID.small.bam -o $workDir/$bamID.small.sorted.bam`;
`samtools index $workDir/$bamID.small.sorted.bam`;
if($options{c}) {
  `samtools depth $workDir/$bamID.small.sorted.bam >$workDir/$bamID.depth`;
}
else {
  `samtools depth -aa -b $workDir/$bamID.bed $workDir/$bamID.small.sorted.bam >$workDir/$bamID.depth`;
}

if(length $prefix) {
  `$prog $workDir $bamID $workDir/$bamID.depth $mt $prefix`;
}
else {
  `$prog $workDir $bamID $workDir/$bamID.depth $mt`;
}

`rm -f $workDir/$bamID.small.bam $workDir/$bamID.small.sorted.bam*`;
if(!$options{i}) {
  `rm -f $workDir/$bamID.bed`;
}
if(!$options{d}) {
  `rm -f $workDir/$bamID.depth`;
}

sub make_bed{
  open(OUT,">$workDir/$bamID.bed") || die "Unable to create $workDir/$bamID.bed\n";
  print OUT "$mt\t1\t16569\n";

  my $rn;
  my $rn1;
  my $region;
  my $str = "$mt:1-16569";

  my @random_set;
  my $range = int ($TotalChromLength/$region_size);
  my %seen;

  for(my $k=0;$k<$N;$k++) {
    my $candidate = int rand($range);
    redo if $seen{$candidate}++;
    push @random_set, $candidate;
  }
  my @sorted_set = sort { $a <=> $b } @random_set;

  my $num_of_regions = 0;
  for(my $k=0;$k<$N;$k++) {
    my $genome_pos = $sorted_set[$k] * $region_size; 
    my $sum_chr_size = 0;
    for(my $i=1;$i<=22;$i++) {
       next if($chr && ($i!=$chr && ("$prefix$i" ne $chr)));
       $rn = $genome_pos - $sum_chr_size;
       $rn1 = $rn+$region_size;
       if($rn>=0 && $rn1<$chr_sizes[$i-1]) {
          print OUT "$prefix$i\t$rn\t$rn1\n";
          $region = "$prefix$i:$rn-$rn1";
          $str = join(' ',$str,$region);
          last;
       }
       $sum_chr_size += $chr_sizes[$i-1];
    }

  }

  close(OUT);

  return $str;
}

sub make_str{
  my ($input_bed) = @_;
  open(IN,"$input_bed") || die "Unable to open $input_bed\n";

  my $str = "$mt:1-16569";
  my $region;
  while(my $line = <IN>)
  {
     chop($line);
     my @b = split('\t', $line);
     next if($b[0] eq $mt || $b[0] eq "MT");
     if($b[0]>=1 && $b[0]<=22) {
       $region = "$prefix$b[0]:$b[1]-$b[2]";
     }
     else {
       $region = "$b[0]:$b[1]-$b[2]";
     }
     $str = join(' ',$str,$region);
  }
  close(IN);

  return $str;
}

#convert the default.bed to be consistent with bam's reference names
sub convert_bed{
  my ($input_bed) = @_;
  open(IN,"$input_bed") || die "Unable to open $input_bed\n";
  open(OUT,">$workDir/$bamID.bed") || die "Unable to create $workDir/$bamID.bed\n";
  print OUT "$mt\t1\t16569\n";

  my $region;
  while(my $line = <IN>)
  {
     chop($line);
     my @b = split('\t', $line);
     next if($b[0] eq "MT");
     if($b[0]>=1 && $b[0]<=22) {
       print OUT "$prefix$b[0]\t$b[1]\t$b[2]\n";
     }
     else {
       print OUT "$b[0]\tb[1]\t$b[2]\n";
     }
  }
  close(IN);
  close(OUT);
}

