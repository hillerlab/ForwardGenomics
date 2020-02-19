#!/usr/bin/perl

# Michael Hiller
# convert a fasta file to a fasta file where all sequences are written in one line
# the parameter -BigProblemChar and -ReplaceChar would replace every BigProblemChar in the sequence by ReplaceChar (was used in the %id pipeline)
# If these parameters are not given, your sequence will be unchanged

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my $BigProblemChar = "";
my $ReplaceChar = "";
GetOptions ("BigProblemChar=s" => \$BigProblemChar, "ReplaceChar=s" => \$ReplaceChar);


$first = 0;
open(file1, $ARGV[0]);
while ($line1 = <file1>) {
	if ($line1 =~/\>.*/) {
		if ($first == 1) {
			print "\n";
		}
		print $line1;
		$first = 1;
	}else{
		chomp $line1;
		if ($BigProblemChar ne "") {
			$line1 =~ s/$BigProblemChar/$ReplaceChar/g;
		}
		print $line1;
	}
}
print("\n");
close(file1);

