#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;

#Command Line Options
my ($inputchrname,$flank,$conversion,$instart,$inend);
my %schema;
my $debug;
my $outdir="";
$flank=0;
my $version = "1.1" ;
my $mean_frag_length=300;
my $minsvsize = 10000;
my $x=4;
my $automean;
&GetOptions(
	'chr=s'=>\$inputchrname,
	'maq|list=s'=>\$conversion,
	'debug!'=>\$debug,
	'start=n'=>\$instart,
	'flank=n'=>\$flank,
	'mean=n'=>\$mean_frag_length,
	'end=n'=>\$inend,
	'minsv=n'=>\$minsvsize,
	'o|dir=s'=>\$outdir,
	'automean!'=>\$automean,
	'x=n'=>\$x,
	'help!'=>\$help
);


my $recordsmatched=0;
if (scalar @ARGV ==0 || $help) {
	system("perldoc $0");
	print "Version:$version\n";
	exit 0;
}

=head1 Synopsis

	Convert Novoalign PE format to BED format for UCSC

	Prints all results to various categories in BED format. If no chr:star-end specified then
	prints everything.
	
	-Multiple novoalign files on command line supported
	-.gz novoalign files are supported.
	- Works best with reads mapped from SEQ/PRB, FASTQ and QSEQ.txt formats


=head2 Command Line Arguments

	An example usage is shown below:

	novope2bed.pl  10-8B.500K.novo  -maq maq.list  -chr chr4 -start  123080000  -end  123820000

	[Fragment Length Settings]
	-mean Set the mean fragment length. Recommended to set this value at (mean + 4 SD). Integer Default = 300.
	-automean Read mean fragment length information from the tail end of the novope output file. Works for .gz and unzipped files. If set will override -mean flag.
		Will add x SDs to the mean and use this as the filter. Boolean switch.		
	-x Set the number of SDs for automean calculation. Integer. Default=4
	
	[Filter Settings]
	-chr chr name. string
	-start start position. integer
	-end end position. integer.
	-flank. Add flanking bp to start and end positions	

	[Preprocessing Settings]	
	-maq maq header conversion file (same as novo2maq). string

	[Output Settings]
	-dir location of output BED files. Defaults to <file.BED> directory. 
	-help  Print this help message

=head2 Authors

Colin Hercus & Zayed Albertyn
Novocraft Technologies & Center for Comparative Genomics  2009

=cut


if ($inputchrname && $instart && $inend) {
	print STDERR  "OPTIONS: $inputchrname $instart-$inend ; Flank= $flank\n";
	$instart-=$flank;
	$inend+=$flank;
	print STDERR  "Seeking : $inputchrname $instart-$inend\n";
}

printf STDERR  "MSG: Input\t%s Files\n\n\n",scalar @ARGV;

foreach $file (@ARGV) {
	next unless -e $file;
	unless ($outdir) {
		$outdir="$file";
	}
	if ($automean) {
		$mean_frag_length = &readMeanFromFile($file,$x);
	}
	processFile($file,$outdir);

}

sub processFile {

	my $file= shift;
	my $outdir = shift;
	&init_color_schema($file,$outdir);
	if ($conversion) {
		my $sedfile = &convert_list($conversion);
		if ($file =~ /\.gz/) {
			open(IN,"zcat $file | sed -f $sedfile |") or die "$!";

		}else {
			open(IN,"sed -f $sedfile $file |") or die "$!";
		}	
	}else {
		if ($file =~ /\.gz/) {
			open(IN,"zcat $file |") or die "$!";
		}else {
			open(IN,$file) or die "$!";
		}
	}

	my $filtered=0;
	my $pairs=0;
	print "MSG: Reading File $file\n";
	while(<IN>) {
		next if /^#/;
		s/>chr/chr/g;
		$pairs++;
		my @F=split(/\s+/,$_);
		#Look at left side only, discard matches without mate aligned
		if ($F[4] eq "U" && $F[10] eq "." && $F[11] ne "." && $F[9] eq "F") {
			my $strand1=$F[9];
			my $strand2=$F[12];
			$strand1 =~ tr/[RF]/[\-\+]/;
			$strand2 =~ tr/[RF]/[\-\+]/;
			my ($start,$end) = ($F[8],$F[11]);
			my $chrname=$F[7];
			$chrname=~s/>//;
			my $name = $F[0];
			$name=~s/\/\d+//;
			my $score=$F[6];
			my $blocksize1=length($F[2]);
		
			my $start1=$start;
			my $end1 = $start + $blocksize1;
		
			my $start2=$end;
			my $end2 = $end +  $blocksize1;

			#filters
			unless ($chrname =~ /^[Cc]hr/) {
				$filtered++;
				next;
			}
			
			
			if($start1 < $start2) {
				$start = $start1;
				$end = $end2;
			} else {
				$start = $start2;
				$end = $end1;
			}
			my $offset = $end - $start - $blocksize1;
			$colour = "0,0,0";
			$colour = "0,255,0" if $offset > $mean_frag_length;
			$colour = "255,0,0" if $strand1 eq $strand2;
			$fhstring = $schema{$colour};
			
			if ($inputchrname && $instart && $inend) {	
			#HERE is where we filter by position of chr and coords supplied from cmd line
			#print "$chrname eq $inputchrname && $start >= $instart && $end <= $inend\n";
				if ($chrname eq $inputchrname && $start >= $instart && $end <= $inend) {
					$recordsmatched++;
					print "MATCHED: $chrname $inputchrname $start<=>$instart $end<=>$inend\n" if $debug;
				}else {
					next;

				}
			}

			if ($strand1 ne $strand2 && (($strand1 eq "+" && $start2 < $start1) || ($strand1 eq "-" && $start1 < $start2)) ) {
				$colour = "255,128,200";
				$fhstring = $schema{$colour};
				#$start -= $blocksize;
				#$end += $blocksize;
				#$offset = $end - $start - $blocksize1;
			}
			next if $strand1 eq $strand2 && $F[1] eq "R";
			$len = $end - $start;
			if($offset <= $blocksize1) {
				$colour = "255,0,255";
				$fhstring = $schema{$colour};
				print  $fhstring "$chrname\t$start\t$end\t$name\t$score\t$strand1\t$start\t$end\t$colour\t1\t$len\t0\n";
				next;
			}
		if ($start1 < $start2 ) {
		       if($offset > $blocksize1 && $offset < $minsvsize) {
					$fhstring = $schema{$colour};
					print $fhstring "$chrname\t$start\t$end\t$name\t$score\t$strand1\t$start\t$end\t$colour\t2\t$blocksize1,$blocksize1\t0,$offset\n";
				}
			}else {
				if($offset > $blocksize1 && $offset < $minsvsize) {
					$fhstring = $schema{$colour};
					print $fhstring "$chrname\t$start\t$end\t$name\t$score\t$strand1\t$start\t$end\t$colour\t2\t$blocksize1,$blocksize1\t0,$offset\n";
				}
			}
			next;
		}
	     if($F[4] eq "U" && ($F[10] ne "." || $F[11] eq ".")) {
			$strand1=$F[9];
			$strand1 =~ tr/[RF]/[\-\+]/;
			$start = $F[8];
			$chrname=$F[7];
			$chrname=~s/>//;
			$name = $F[0];
			$name=~s/\/\d+//;
			$score=$F[6];
			my $blocksize=length($F[2]);
		
			$end = $start + $blocksize;
		
			$colour = "0,0,255";
			$colour = "0,255,255" if $F[11] eq ".";
		        if ($inputchrname && $instart && $inend) {      
				#HERE is where we filter by position of chr and coords supplied from cmd line
				if ($chrname eq $inputchrname && $start >= $instart && $end <= $inend) {
					$recordsmatched++;
					print "MATCHED: $chrname $inputchrname $start<=>$instart $end<=>$inend\n" if $debug;
				}else {
					next;
				}
			}
			$fhstring = $schema{$colour};
			print $fhstring "$chrname\t$start\t$end\t$name\t$score\t$strand1\t$start\t$end\t$colour\t1\t$blocksize\t0\n";
		}

	}

	print "SUMMARY:\tNo. mate-pairs read: $pairs\n";
	print "SUMMARY:\tFiltered out non-chr names: $filtered\n";
	print "SUMMARY:\tMatched input: $recordsmatched\n";

	close IN;
	print "SUMMARY:Record Counts\n";
	#Print the line count summaries
	foreach $outfile (values %schema) {
		my $count= `wc -l $outfile.out.bed`;
		close $outfile;
		my($number,$filename) =split(/\s+/,$count);
		printf "$filename\t%s\n",$number-1;
		
		if ($inputchrname && $instart && $inend) {
			my $outdir ="$inputchrname.$instart-$inend";
			mkdir($outdir) unless -e $outdir;
			system("mv $outfile.out.bed $outdir/ ");
		}
	}

}

#this is where i control colour schema and open a separate file for each 
sub init_color_schema {
		my $file =shift;
		my $outdir=shift;
		
		$outdir=~s/\///g;
		$outdir.=".BED";
		if (! -e $outdir ) {
			mkdir $outdir;
		}
	
		$file=~s/\.novo|\.novop//g;
		$file="$outdir/$file";
		$schema{"0,0,0"}="$file.normal";
		$schema{"0,255,0"}="$file.long";	
		$schema{"255,0,0"}="$file.inverted";	
		$schema{"255,0,255"}="$file.overlap";	
		$schema{"0,255,255"}="$file.singleton";	
		$schema{"0,0,255"}="$file.chimera";	
		$schema{"255,128,200"}="$file.negative";
		my $name;
		my $colour;
		foreach $colour ( keys %schema) {
			$name=$schema{$colour};
			open($name,"+>$name.out.bed");
			print $name "track name=\"$name\" description=\"$name \" color=$colour\n";
			#print STDERR "MSG: opening $name.out.bed\n";
		}

}


sub convert_list {
	my $file =shift;
	open(IN,$file) or die "$!";
	my $sedfile ="$file.sed";
	open(SED,"+>$sedfile");
	while(<IN>) {
		chomp;
		s/>//g;
		my($name,$newname) = split;
		print SED "s/$name/$newname/;\n";
	}
	close IN;
	close SED;
	print STDERR  "MSG:Successfully made $sedfile for sed\n";
	return $sedfile;
}


sub readMeanFromFile {
	my $file=shift;
	my $x=shift;
	my $mean=0;
	my $sd=0;
	my $pipe="";
	if ($file =~ /\.gz$/) {
		$pipe="gunzip -c $file | tail -n 100 |";
	}else {
		$pipe = "tail -n 100 $file |";
	}
	open(PI,$pipe) or die "$!";
	while(<PI>) {
		if (/^# *Mean *([0-9]+), *Std Dev *([0-9.]+)/) {
			$mean=$1;
			$sd=$2;
		}
	}
	close PI; 
	print STDERR "$file: Auto Mean filter: ($mean + $x x $sd)\n";
	return $mean + ($x*$sd);

}


