#!/usr/bin/perl

# Michael Hiller, 2015
# simple script to convert a tree.mod into a newick file. 
# If the input is already a newick file, it will output the same tree (the tree must be given on a single line)

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

$| = 1;		# == fflush(stdout)
# options
my $verbose = 0;

sub usage {
	die "usage: $0 in.treeFile out.treeFile\n\tin.treeFile is either a tree in phastCons model format or newick.\n\tout.treeFile will be in newick format\n\tNOTE: If you input a newick tree already, the tree must be given on a single line\n";
};
GetOptions ("v|verbose"  => \$verbose) || usage();
usage() if ($#ARGV < 1);


open(FILEOUT, ">$ARGV[1]") || die "ERROR: cannot write to $ARGV[1]\n";

open(FILEIN, "$ARGV[0]") || die "ERROR: cannot open to $ARGV[0]\n";
my @content = <FILEIN>;
close FILEIN;


my $success = 0;
for (my $i=0; $i<=$#content; $i++){
	my $line = $content[$i];
	if ($line =~ /^TREE: (.*)/) {
		print "found TREE line: $line\n" if ($verbose);
		checkAndOutput($1);
	}elsif ($line =~ /^\(/) {
		print "found newick input: $line\n" if ($verbose);
		checkAndOutput($line);
	}else{
		print "skip line: $line\n" if ($verbose);
	}
}
close FILEOUT;

die "ERROR: could not find a newick tree in $ARGV[0]\n" if ($success == 0);
exit 0;



############################
# check number of ( and ) and if they match, output the tree
############################
sub checkAndOutput {
	my $line = shift;
	
	my $count1 = ($line =~ s/\(/\(/g);
	my $count2 = ($line =~ s/\)/\)/g);
	die "ERROR: found the newick string but the number of opening and closing parenthesis don't match ($count1 != $count2) : $line\n" if ($count1 != $count2);
	print FILEOUT "$line";
	print "output correct newick format: $line\n" if ($verbose);
	$success = 1;
}
