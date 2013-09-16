#!/usr/bin/perl -w

use File::Basename;
use Data::Dumper;
use Getopt::Long;

#Copyright Zayed Albertyn

my ($gzip,$file,$k,$step,$norun,$paf,$build,$help,$format,$opts,$stdout,$db);


&GetOptions(
	'formatdb|build|index!' => \$build,
	'norun!'=>\$norun,
	'k=n'=>\$k,
	's=n'=>\$step,
	'help!'=>\$help,
	'zip|compress!'=>\$gzip,
	'pafeval|paf!'=>\$paf,
	'db=s'=>\$db,
	'opts=s'=>\$opts,
	'format=s'=>\$format,
	'stdout!'=> \$stdout
);

unless ($opts) {
	$opts="";
}


if (scalar @ARGV ==0 || $help) {
	system("perldoc $0");
	exit 1;
}

unless ($format) {
	$format="Eland";
}


if ($build && $db) {
	&Build;
	print STDERR "$db built\nExiting\n";
	exit 0;
}


unless (-e $db) {
	system("perldoc $0");
	exit;
}
my @SUFFIXES = qw(.fastq .fasta .fa .fna .txt .fa);


&Run;





#run the program on each file
sub Run {
	my $ext;
	if ($opts) {
		my $optstring=$opts;
		$optstring=~s/\-|\s+//g;
		$ext="novo_$optstring";
	}else {
		$ext="novo";	
	}

	my $tag="-q";
	foreach (@ARGV){
		my $path=$_;
		 my($filename, $directories, $suffix) = fileparse($path, @SUFFIXES);
	
		if (/_sequence.txt$|fastq$/) {
			$tag="-f"
		}
		elsif (/\.fa$|\.fasta$|\.fna/) {
			$tag="-f";
		}
		elsif (/prb.txt/) {
			$tag="-q";
		}
		my $outfile="$filename.$ext";	
		my $cmd = "novoalign $tag $_ -d $db -o $format  $opts";
		print "#$cmd > $outfile\n";
		
		unless ($norun) {
			system("gtime -v -o $outfile.time $cmd  > $outfile");

			if ($paf && -e $outfile && $format eq "Eland") {
					print "# paf_utils.pl pafeval -p $outfile.paf\n";
					system("paf_utils.pl novo2paf $outfile > $outfile.paf");
					system("paf_utils.pl pafeval $outfile.paf ");
			}

		}
		if ($gzip) {
			system("gzip -f  $outfile $filename.paf");
		}

			
	}

}

#Build the database
sub Build {
	system("novoindex $db.nv $db");0
}


sub Docs {

=head1 SYNOPSIS

	A wrapper for running novoalign

=head2 Command Line Options

		-db Database name. can also be fasta file that will be indexed
		-build index the database
		-norun dont run anything, just print the commands.Useful for checking large jobs 
		-format output format. Choices are Eland (default), Novo or Blast
		-opts String of options to supply to novoalign. run novoalign without args to see these
		
		-pafeval Convert Eland-formatted output to PAF and run pafeval
			 Note that this option should only be used on simulated datasets where the read naming convention
			 contains the original sequence location
		
		-del  Delete the novo output file 

=head2 Author: Zayed Albertyn


=cut





}

