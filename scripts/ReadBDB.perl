#!/sw/bin/perl

# Michael Hiller, 2015
# just extracts the value for a key from a given BDB
# if no key is given, print entire content 
# option to output every key into a separate file $outputDir/$key.$suffix

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BerkeleyDB;

my $brief = 0;				# no additional key output and separating newlines; just the values
my $outputDir = "";
my $suffix = "";
sub usage {
	die "usage: $0 INPUT.bdb [key] [-brief] [-outputDir string -suffix string]
# Without any parameters, it outputs the entire BDB file in the format:
value for El10993:
VALUE ....
value for El10994:
VALUE ....
\
# If you specify a key (El10993), it will only output the value of this key
\
# If -brief is given, the output will not contain the \"value for El10994:\" lines
\
# If -outputDir DIR -suffix SUFFIX is given, it will try to create this output dir and then create a DIR/ELEMENT.SUFFIX file with the value for every key. 
# If you additionally specify a single key, you will only get the value of this key in DIR/ELEMENT.SUFFIX.
# You must set both -outputDir and -suffix. 
# Example: $0 Alignment.bdb -outputDir in -suffix fa
\n"
};
GetOptions ("brief" => \$brief, "outputDir=s" => \$outputDir, "suffix=s" => \$suffix) || usage();
usage() if ($#ARGV < 0);

die "ERROR: -outputDir is given. You must also specify the suffix of the output file with -suffix SUFFIX\n" if ($outputDir ne "" && $suffix eq "");
die "ERROR: -suffix is given but outputDir is not set. If you want to get the values in one file per key, set -outputDir DIR\n" if ($outputDir eq "" && $suffix ne "");


# create output dir if necessary
if ($outputDir ne "") {
	my $call = "set -o pipefail; mkdir -p $outputDir";
	`$call`;
	die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
}


# open BDB
my %h;
tie %h, "BerkeleyDB::Btree",	-Filename => "$ARGV[0]", -Flags => DB_RDONLY || die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";


# read from the BDB
# output to stdout
if ($outputDir eq "") {
	if ($#ARGV > 0) {
		print "value for $ARGV[1]\n" if (! $brief);
		if (exists $h{$ARGV[1]}) {
			print "$h{$ARGV[1]}";	
		}else{
			print "ERROR: KEY $ARGV[1] does not exist\n";
		}
	}else{
		foreach my $key (keys %h) { 
			print "value for $key:\n" if (! $brief);
			if (exists $h{$key}) {
				print "$h{$key}";	
			}else{
				print "ERROR: KEY $key does not exist\n";
			}
			print "\n" if (! $brief);
		}
	}
# output to directory
}else{
	if ($#ARGV > 0) {
		if (! exists $h{$ARGV[1]}) {
			print "ERROR: KEY $ARGV[1] does not exist\n";
		}else{
			my $file = "$outputDir/$ARGV[1].$suffix";
			open FILEOUT, ">$file" || die "ERROR: cannot write to $file\n";
			print FILEOUT "$h{$ARGV[1]}";	
			close FILEOUT;
		}
	}else{
		foreach my $key (keys %h) { 
			if (exists $h{$key}) {
				my $file = "$outputDir/$key.$suffix";
				open FILEOUT, ">$file" || die "ERROR: cannot write to $file\n";
				print FILEOUT "$h{$key}";	
				close FILEOUT;
			}else{
				print "ERROR: KEY $key does not exist\n";
			}
		}
	}
}

# close BDB
untie %h ;

