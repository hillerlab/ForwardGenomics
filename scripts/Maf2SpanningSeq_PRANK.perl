#!/usr/bin/perl

# Michael Hiller, 2015
# given a maf file listing a single region (cons el), read all maf blocks and determine the sequence coordinates that cover maf block 1-N
# NOTE: These coords include large insertions and unaligning regions in the middle, but not unaligning regions (e-line blocks) at both ends
# Only e-line blocks will be included that are spanned by s-line blocks.
# Then extract these coords from the 2bit Files. 
# Optionally, 
#     - run prank
#     - output alignments of species and named ancestors 
#     - add the species (and internal nodes in case some are sister species) that have a complete deletion 
# Prank reproducibly fails if the base composition of the input is highly biased (typically AG rich). We catch these errors and re-run prank setting the base composition to 60%AT. 
# The script cleans up temp files before error-exit. 

use strict;
use warnings;
use Env;
use lib "./";

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
#use Globals;
use MyKentFunctions;
use MyBDB;

my $twoBitPath = "";

# options
# Pls see the usage() for an explanation of each parameter
my $verbose = 0;
my $TwoBitSuffix = ".2bit";
my $BDBFile = "";
my $noMissingSeq = 0;
my $treeFile = "";
my $runPrankFlag = 0;
my $keepTemporaryFiles = 0;
my $maxAlignmentLength = 10000000000000;		# do not run prank if the total number of bp in the ali exceed this number. Purpose: dont start long jobs until the user pushes them to the long queue

sub usage {
	die "usage: $0 Maf_or_BDBFile ElementID -twoBitPath  -v|-verbose  [PARAMETERS]
\
$0 will read the maf from either a BDB file or a maf file (must have the suffix .maf) produced by mafExtract
It gets the sequences spanned by the first and last s-line block for each species and outputs them.
The -twoBitPath parameter must point to the directory that contains the species/species.2bit files. 
Optionally, you can directly run prank and $0 outputs the alignment and the named ancestors. 	
\
general parameters:
	-TwoBitSuffix string [default .2bit]      suffix of the twoBitPath/species/species.$TwoBitSuffix 2bit file. Use .quality.2bit to use the quality-masked files
	-noMissingSeq                             flag: if set, we will exclude species that do not have s or e lines for each given maf block (the missing blocks must be assembly gaps)
\
parameters for running prank:
	-runPrank                                 flag: if set, run prank and output the alignment + the ancestors
	-treeFile string                          file with the tree in either phastCons mod or in newick format. MANDATORY if -runPrank. 
	-BDBFile string                           name of the BDBFile to output the result. Only necessary if -runPrank is given.
	-keepTemporaryFiles                       flag: keep the prank files.
	-maxAlignmentLength int                   error-exit if the number of bp in the ali exceed this threshold (default $maxAlignmentLength). Use it to filter out long-running jobs for the long queue.\n"
};
GetOptions ("twoBitPath=s" => \$twoBitPath, "TwoBitSuffix=s" => \$TwoBitSuffix, "v|verbose"  => \$verbose,	"runPrank" => \$runPrankFlag, 
	"BDBFile=s" => \$BDBFile, "treeFile=s" => \$treeFile, "maxAlignmentLength=i" => \$maxAlignmentLength,
	"keepTemporaryFiles" => \$keepTemporaryFiles, "noMissingSeq" => \$noMissingSeq) || usage();
usage() if ($#ARGV < 1);

die "ERROR: You must set the -twoBitPath parameter\n" if ($twoBitPath eq "");

if ($verbose) {
	$| = 1;		# == fflush(stdout)
}

# if set to 1, convert all lower case letters in the 2bit file to upper case. Use this is to unmask repeat-masked seq from a regular 2bit file. 
# if set to 2, convert lower case to N. We do this if the TwoBitSuffix contains "quality" indicating our quality-masked 2bit files where low quality = lower case. 
my $unmaskFlag = 1;
if ($TwoBitSuffix =~ /quality/i) {
	$unmaskFlag = 2;
	print "TwoBitSuffix contains the word \"quality\" --> will replace lower case characters from the 2bit file by N's\n";
}

my $ElementID = $ARGV[1];

# global hashes and variables
my %aliSpeciesAndAnc;						# hash that will have the final aligned sequences of all species (without missing data) and all ancestors
my %skipSpecies;								# hash of species that are given in the maf but where some maf blocks have missing sequence --> we will skip them
my %completeDeletionSpecies;				# hash of species where the entire element is deleted or unaligning (only e lines in the maf) --> we will add the "----..." sequence later
my %speciesWithoutMissingData;			# hash of species where all maf blocks are given (sequence and or deletion/unaligning) --> all these species will appear in the final output
my %species2Chrom;							# hashes used to parse the given maf
my %species2Start;
my %species2End;
my %species2StartSblock;
my %species2EndSblock;
my %species2Strand;
my %species2srcSize;
my %species2NumBlocksWithData;			# for how many maf blocks does a species have e or s lines
my $totalBlockNum = 0;						# total number of maf blocks
my %species2NumSBlocks;

# if $runPrankFlag test if BDBfile and the tree is given
if ($runPrankFlag) {
	die "ERROR: -runPrank is given, but parameter -treeFile filename is not given \n" if ($treeFile eq "");
	die "ERROR: -runPrank is given, but the tree file does not exist (either -treeFile is not given or the default $treeFile does not exist)\n" if (! -e $treeFile);
	die "ERROR: -runPrank is given, but flag -BDBFile is not given\n" if ($BDBFile eq "");
}
doit($ARGV[0]);		# just that file


##################################################################
# read maf from file or extrac it from BDB, get the spanning seq and optionally run prank 
##################################################################
sub doit {
	my $file = shift;
	
	# read the maf from a file or BDB
	my @maf;
	if ($file =~ /\.maf$/) {
		open(file1, $file) || die "ERROR: cannot open $file\n";
		@maf = <file1>;
		close file1;
		chomp(@maf);
	}else{
		my $mafFromBDB = readBDB($file, $ElementID); 
		@maf = split(/\n/, $mafFromBDB);
	}
	
	if ($verbose) {
		$" = "\n";
		print "maf read from $file:\n@maf\n\n";
	}

	# parse the maf. First count the number of mafblocks and get the reference species (species in the first s-line)
	my $referenceSpecies = "";
	for (my $i=0; $i<=$#maf; $i++){
		my $l = $maf[$i];
		if ($l =~ /^a score/) {
			$totalBlockNum ++;
		}elsif ($referenceSpecies eq "") {
			if ($l =~ /^s .*/) {
				my ($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($l);
				$referenceSpecies = $species;
			}
		}
	}
	print "--> reference species $referenceSpecies   #mafblocks: $totalBlockNum\n" if ($verbose);

	# now read the s/e/i lines 
	my $curBlockNum = 0;
	for (my $i=0; $i<=$#maf; $i++){
		my $l = $maf[$i];
		
		if ($l =~ /^a score/) {
			$curBlockNum ++;

		}elsif ($l =~ /^[se] .*/) {		# parse s or e lines

			my ($species, $chr, $start, $size, $strand, $srcSize, $seq, $status);
			if ($l =~ /^s .*/) {
				($species, $chr, $start, $size, $strand, $srcSize, $seq) = getSLineData($l);
				$species2NumSBlocks{$species} ++;

				# read the following i-line
				# only the reference species have no i-lines
				if ($species ne $referenceSpecies) {
					$i++;
					my $iline = $maf[$i];
					if ($iline =~ /^i .*/) {
						my ($species, $chr, $statusUp, $countUp, $statusDown, $countDown) = getILineData($iline);
						# If there is not a C or I status upstream, then the sequence is not continous (meaning a discontinuity) --> skip
						# not for the first block
						if ($curBlockNum > 1) {
							if ($statusUp ne "C" && $statusUp ne "I") {
								print STDERR "$ElementID $species has a discontinuity: status upstream is $statusUp for block $curBlockNum ($iline)\n--> will skip $species\n";
								$skipSpecies{$species} = 1;
							}
						}
						if ($curBlockNum < $totalBlockNum) {
							if ($statusDown ne "C" && $statusDown ne "I") {
								print STDERR "$ElementID $species has a discontinuity: status downstream is $statusDown for block $curBlockNum ($iline)\n--> will skip $species\n";
								$skipSpecies{$species} = 1;
							}
						}
					}else{
						die "ERROR: there is no i-line that follows this s-line in the maf:\n\t$l\n\t$iline\n";
					}
				}
				
			}elsif ($l =~ /^e .*/) {
				($species, $chr, $start, $size, $strand, $srcSize, $status) = getELineData($l);
				$species2NumSBlocks{$species} = 0 if (! exists $species2NumSBlocks{$species});
			}
			
			my $end = $start+$size;
			# check if the chrom, strand is the same and if the order of the blocks is consistent (on both + and - strand, the coords for block start and end are increasing)
			if (exists $species2Chrom{$species} && $species2Chrom{$species} ne $chr) {
				print STDERR "$species has $species2Chrom{$species} as the chrom from previous s-lines, which differs from $chr in ($l)\n--> will skip $species because maf blocks come from different chroms\n";
				$skipSpecies{$species} = 1;
			}
			if (exists $species2Strand{$species} && $species2Strand{$species} ne $strand) {
				print STDERR "$species has $species2Strand{$species} as the strand from previous s-lines, which differs from $strand in ($l)\n--> will skip $species because maf blocks come from different strands\n";
				$skipSpecies{$species} = 1;
			}
			if (exists $species2Start{$species} && $species2Start{$species} > $start) {
				print STDERR "$species has $species2Start{$species} as the start from previous s-lines, which is >= start $start in ($l)\n--> will skip $species because maf block starts seem to be mis-ordered\n";
				$skipSpecies{$species} = 1;
			}
			if (exists $species2End{$species} && $species2End{$species} > $end) {
				print STDERR "$species has $species2End{$species} as the end from previous s-lines, which is >= end $end in ($l)\n--> will skip $species because maf block ends seem to be mis-ordered\n";
				$skipSpecies{$species} = 1;
			}

			# store the start/end of the first/last s-line block, as we will only extract sequence between s-line blocks (flanking e line blocks will be ignored)
 			if ($l =~ /^s .*/) {
				# start --> only the first block
				$species2StartSblock{$species} = $start if (! exists $species2StartSblock{$species});
				# always update the end
				$species2EndSblock{$species} = $end;
			}

			# store the data
			$species2Chrom{$species} = $chr if (! exists $species2Chrom{$species});
			$species2Strand{$species} = $strand if (! exists $species2Strand{$species});
			$species2Start{$species} = $start if (! exists $species2Start{$species});
			$species2End{$species} = $end;		# ends are always updated 
			$species2srcSize{"$species.$chr"} = $srcSize;
			$species2NumBlocksWithData{$species} = 0 if (! exists $species2NumBlocksWithData{$species});  # init
			$species2NumBlocksWithData{$species} ++; 

		}
	}

	die "ERROR: no maf block given in $file for $ElementID\n" if ($totalBlockNum == 0);

	# if -noMissingSeq is set, check which species don't have data for each maf block and add them to the skipSpecies hash
	if ($noMissingSeq == 1) {
		foreach my $species (keys %species2Chrom) {
			next if (exists $skipSpecies{$species});
			if ($species2NumBlocksWithData{$species} != $totalBlockNum) {
				print "\t$species is skipped because it has only data for $species2NumBlocksWithData{$species} of the total of $totalBlockNum maf blocks (-noMissingSeq is set)\n";
				$skipSpecies{$species} = 1;
			}
		}
	}
	
	if ($verbose) {
		print "-- spanning coordinates after reading $totalBlockNum blocks\n";
		outputSpanningCoords();
	}

	# now run prank or just output the seqs
	if ($runPrankFlag) {
		# run prank and output ali + ancestors
		runPrank();
	}else{
		# output the seqs
		outputSpanningSeq("");
	}
	
	return 0;
}


#############################################################
# convert the given - to + strand coords, using srcSize
#############################################################
sub MinusToPlusStrandCoords {
	my ($start, $end, $srcSize) = @_;
	
	my $startPlus = $srcSize-$end;
	my $endPlus = $srcSize-$start;
	
	# sanity check
	die "ERROR in MinusToPlusStrandCoords: cannot convert ($start, $end, $srcSize) to +strand coords as this exceeds srcSize (+ coords: $startPlus - $endPlus)\n" if ($startPlus >= $srcSize || $endPlus > $srcSize);
	die "ERROR in MinusToPlusStrandCoords: cannot convert ($start, $end, $srcSize) to +strand coords as results in <0 coords (+ coords: $startPlus - $endPlus)\n" if ($startPlus < 0 || $endPlus < 0);
		
	return ($startPlus, $endPlus);
}


#######################################################################
# only if -verbose is set : small function to output the spanning coords
#######################################################################
sub outputSpanningCoords {
	my $skippedList = "";
	foreach my $species (keys %species2Chrom) {
		if (exists $skipSpecies{$species}) {
			$skippedList .= "\t$species is skipped because of inconsistencies (different chrom, strand, coordinate order)\n";
			next; 
		}

		# fake init for those species who have only e-line blocks. Just to avoid uninitialized hash values below.
		if ($species2NumSBlocks{$species} == 0) {
			$species2EndSblock{$species} = $species2End{$species};
			$species2StartSblock{$species} = $species2Start{$species};
		}

		my $span = $species2EndSblock{$species}-$species2StartSblock{$species};
		printf "\t%12s\t", $species;
		print "$species2Chrom{$species}:$species2StartSblock{$species}-$species2EndSblock{$species} $species2Strand{$species}\t\tsrcSize ", $species2srcSize{"$species.$species2Chrom{$species}"}, "\t\tspan $span\t\t";
		if ($species2Strand{$species} eq "-") {
			my ($startPlus, $endPlus) = MinusToPlusStrandCoords($species2StartSblock{$species},$species2EndSblock{$species},$species2srcSize{"$species.$species2Chrom{$species}"});
			print "+coords: $species2Chrom{$species}:$startPlus-$endPlus";
		}
		print "  !!!! ONLY E-line blocks given --> seq will not be in the output !!!" if ($species2NumSBlocks{$species} == 0);
		print "\n";
	}
	print "$skippedList\n";
}

#######################################################################
# small function to get the spanning seqs using our twoBitToFa wrapper function
# If a non-empty filename is given, the output will happen to this file
#######################################################################
sub outputSpanningSeq {
	my $filename = shift;
	
	# count the number of species in the fasta file and return it. If there is only one species than there is nothing to align. 
	my $speciesCount = 0;
	
	if ($filename ne "") {
		open FILEOUT, ">$filename" || die "ERROR: cannot write to file $filename\n";
	}
	
	my $totalNumBases = 0;
	foreach my $species (keys %species2Chrom) {
		if (exists $skipSpecies{$species}) {
			next; 
		}
		# add this species to the list of species without missing data (sequence and or complete deletion)
		$speciesWithoutMissingData{$species} = 1;

		# do not output the seq if only e-line blocks are given (this completely unrelated seq hurts the prank alignment). 
		# Also if it is a deletion, then there is no sequence to give to prank
		# We add these species to the hash and add the "------..." sequence later back to the alignment
 		if ($species2NumSBlocks{$species} == 0) {
			$completeDeletionSpecies{$species} = 1;
			next;
		}

		# convert - to + strand 
		my $startPlus = $species2StartSblock{$species};
		my $endPlus = $species2EndSblock{$species};
		if ($species2Strand{$species} eq "-") {
			($startPlus, $endPlus) = MinusToPlusStrandCoords($species2StartSblock{$species},$species2EndSblock{$species},$species2srcSize{"$species.$species2Chrom{$species}"});
		}
		
		# retrieve the seq and print
		my $seq = getSeqTwoBitFile($species2Chrom{$species}, $startPlus, $endPlus, $species2Strand{$species}, "$twoBitPath/$species/$species$TwoBitSuffix", $unmaskFlag);
		print "retrieve ", $endPlus-$startPlus, " bp sequence for $species $species2Chrom{$species}:$startPlus-$endPlus from 2bit file: $twoBitPath/$species/$species$TwoBitSuffix\n";
		if ($filename ne "") {
			print FILEOUT ">$species\n$seq\n";
			$speciesCount++;
		}else{
			print ">$species\n$seq\n";
		}
		$totalNumBases += $endPlus-$startPlus; 
	}
	print "total number of bases to be aligned for $ElementID: $totalNumBases \n";
	
	# close file if necessary
	if ($filename ne "") {
		close FILEOUT;
	}
	
	return ($speciesCount, $totalNumBases);
}


#######################################################################
# remove tmpfiles and die
#######################################################################
sub cleanDie {
	my ($message, $tmpFile) = @_;
	`rm -f ${tmpFile}*` if (! $keepTemporaryFiles);
	die($message);
}

#######################################################################
# run the given prank call. 
# Sometimes prank crashes here with an error "Reading the alignment file failed". 
#    Ari said this comes from the empirically estimated base freqs and indeed we had a distortion of e.g. A and G rich seqs
#    A workaround is to set the basefreq (here to the genomewide average 60% AT) and according to Ari
#    "would that affect the ancestral reconstruction in anyway? No, reconstruction of ancestral character states is done with BppAncestor and Prank settings don't affect that."
#    This turns out to be not true, as Michael found examples where -dnafreqs makes a minor difference. 
#    However, the differences are rare and minor, and we need a method that is stable.  
# Therefore, we run prank and only if it fails with this error, we set the -dnafreqs parameter
#
# The expected $call should not redirect the stdout of prank, as we are doing this in this function explicitely.
#######################################################################
sub runPrankAndCheckForError {
	my ($call, $tmpFile) = @_;
	
	# generate the call with -dnafreqs and redirect stdout in both calls
	my $callWithBaseFreqs = "$call -dnafreqs=30,20,20,30 > ${tmpFile}.prankout";
	$call = "$call > ${tmpFile}.prankout";
	
	print "Running prank: $call   ...\n";
	system("$call");
	if ($? == 0) {
		print "DONE\n";
		return 0;
	}else{
		# see if the output file contains the error message
		my $grepError = "grep \"Reading the alignment file failed. Exiting\" ${tmpFile}.prankout -c";
		my $result = `$grepError`;
		cleanDie("ERROR: $grepError failed\n", $tmpFile) if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
		chomp($result);
		if ($result > 0) {
			print "WARNING:\tPrank failed with error message:  Reading the alignment file failed. Exiting\n\t\t	Now set -dnafreqs and run again\n";
			print "Running prank: $callWithBaseFreqs   ...\n";
			system("$callWithBaseFreqs");
			if ($? == 0) {
				print "DONE\n";
				return 0;
			}else{
				cleanDie("ERROR: $callWithBaseFreqs failed\n", $tmpFile);
			}
 		}else{
			cleanDie("ERROR: $call failed\n", $tmpFile);
		}
	}
	return 0;
}

#######################################################################
# put the seqs into a temp.fa file, run prank, parse the output and return the alignment & the ancestral sequences
#######################################################################
sub runPrank {

	# use /dev/shm for temp files as lustre has issues with creating and deleting many files
# 	my $tmpFile = `mktemp $ENV{'TMPDIR'}/prank.in.XXXXXXXXX.fa`;
	my $tmpFile = `mktemp /dev/shm/prank.in.XXXXXXXXX.fa`;
	chomp($tmpFile);
	print "Using tempfile $tmpFile ...\n" if ($verbose);

	# get the sequences into the tmpFile
	my ($speciesCount, $alignmentLength)  = outputSpanningSeq($tmpFile);
	
	if ($alignmentLength > $maxAlignmentLength) {
		print STDERR "WARNING_maxAlignmentLength:\t$ElementID\tMaximum Alignment Length of $maxAlignmentLength is exceeded (total number of bases: $alignmentLength). Will not start prank. Push $ElementID to the long queue and do not specify -maxAlignmentLength\n";
		secureWriteBDB($BDBFile, $ElementID, "NOTHING_TO_ALIGN\n", 1, ".");
		`rm -f ${tmpFile}*` if (! $keepTemporaryFiles);
		exit(0);
	}
	
	if ($speciesCount < 2) {
		print "only $speciesCount species in the fasta file ($tmpFile) --> nothing to align --> add string NOTHING_TO_ALIGN to $BDBFile\n";
		secureWriteBDB($BDBFile, $ElementID, "NOTHING_TO_ALIGN\n", 1, ".");
		# cleanup, unless $keepTemporaryFiles is set
		`rm -f ${tmpFile}*` if (! $keepTemporaryFiles);
		print "*** ALL DONE ***\n";
		return;
	}

	# convert the tree to newick format (regardless of whether it is mod or nh format)
	my $call = "set -o pipefail; ConvertTreeModToNewick.perl $treeFile ${tmpFile}.tree.nh";
	system("$call") == 0 || cleanDie("ERROR: $call failed\n", $tmpFile);
	if ($verbose) {
		print "\ntreeFile converted to newick: ${tmpFile}.tree.nh\n";
		system "cat ${tmpFile}.tree.nh";
		print "\n\n";
	}

	# this will call prank, output files will be $tmpFile.best.anc.fas etc
	# NOTE: There is a random element in prank of where to place ambiguous indels. Therefore, lets fix the seed so that we always get the same.
	$call = "set -o pipefail; prank -d=$tmpFile -t=${tmpFile}.tree.nh -o=$tmpFile -showtree -showanc -prunetree -once -seed=10";
	runPrankAndCheckForError($call, $tmpFile);

	# now add the species (and internal nodes if necessary) that have complete deletions
	# We do this by adding those complete-deletion species to the prank alignment (only species) and then run prank again to estimate only the ancestors
	my $removeFirstAlignmentColumn = addCompleteDeletions($tmpFile);
	# the result is ${tmpFile}.best.anc.dnd and ${tmpFile}.best.anc.fas  (regardless of whether we ran prank again (because there are complete-deletion species) or not

	# now create a map between the #2# ancestor names and our mm10-rn5 names
	# (((((((mm10:0.08360,rn5:0.08948)#1#:0.22169,speTri2:0.13675)#7#:0.00958,cavPor3:0.23058)#11#:0.02794,oryCun2:0.21242)#14#:0.01413,
	# (((((((mm10:0.0836,rn5:0.08948)mm10-rn5:0.22169,speTri2:0.13675)mm10-speTri2:0.00958,cavPor3:0.23058)mm10-cavPor3:0.02794,oryCun2:0.21242)mm10-oryCun2:0.01413,
	# #1# --> mm10-rn5
	# #7# --> mm10-speTri2 
	# etc
	#
	# first remove all #[0-9]*# from the ancestor tree and then run tree_doctor to name the ancestors
	$call = "set -o pipefail; cat ${tmpFile}.best.anc.dnd | sed 's/#[0-9]*#//g' | tree_doctor /dev/stdin -a -n";
	print "$call\n" if ($verbose);
	my $namedTree = `$call`;
	cleanDie("ERROR: $call failed\n", $tmpFile) if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);	
	chomp($namedTree);
	# now read the prank tree (NOTE: This tree is identical, except for the anc names
	$call = "set -o pipefail; cat ${tmpFile}.best.anc.dnd";
	print "$call\n" if ($verbose);
	my $prankTree = `$call`;
	cleanDie("ERROR: $call failed\n", $tmpFile) if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);	
	chomp($prankTree);
	# This will convert 
	# (((((((mm10:0.08360,rn5:0.08948)#1#:0.22169,speTri2:0.13675)#7#:0.00958,cavPor3:0.23058)#11#:0.02794,oryCun2:0.21242)#14#:0.01413,
	# into 
	# mm10 rn5 #1# speTri2 #7# cavPor3 #11# oryCun2 #14#
	$prankTree =~ s/:[0-9.]*[,;)]/ /g; 
	$prankTree =~ s/\(//g;
	$prankTree =~ s/;//g;
	$namedTree =~ s/:[0-9.]*[,;)]/ /g; 
	$namedTree =~ s/\(//g;
	$namedTree =~ s/;//g;
	print "\tprankTree: $prankTree\n\tnamedTree: $namedTree\n" if ($verbose);
	# now split and create the map
	my @prankNodes = split(/[ ]+/, $prankTree);
	my @namedNodes = split(/[ ]+/, $namedTree);
	cleanDie("ERROR in ancestral name mapping: Number of fields differ between $prankTree and $namedTree\n", $tmpFile) if ($#prankNodes != $#namedNodes);
	my %AncNum2AncName;
	for (my $i=0; $i<=$#prankNodes; $i++) {
		if ($prankNodes[$i] =~ /#\d+#/) { 
			print "\tmap: $prankNodes[$i] --> $namedNodes[$i]\n" if ($verbose);
			$AncNum2AncName{$prankNodes[$i]} = $namedNodes[$i];
		}else{
			cleanDie("ERROR in ancestral name mapping: non-ancestral nodes differ in their name: $prankNodes[$i] vs. $namedNodes[$i]\n", $tmpFile) if ($prankNodes[$i] ne $namedNodes[$i]);
		}
	}
	
	# now get the aligned species and ancestors (after renaming with their proper name) in oneline-fasta format
	my $finalAlignment = "";
	open FILEIN, "${tmpFile}.best.anc.fas" || cleanDie("ERROR: cannot read from ${tmpFile}.best.anc.fas\n", $tmpFile);
	my @ali = <FILEIN>;
	close FILEIN;
	chomp(@ali);
	my $currentName = "";		# proper name of the species or internal node
	for (my $i=0; $i<=$#ali; $i++) {
		# fasta header vs sequence
		if ($ali[$i] =~ /^>(.*)/) {
			my $name = $1;
			# if the name starts with # use the map to output named ancestors
			if ($name =~ /^#/) {
				cleanDie("ERROR in writing final results: Cannot find ancestor $name in the AncNum2AncName hash\n", $tmpFile) if (! exists $AncNum2AncName{$name});
				# always add a newline in front of a new fasta header, unless this is the very first 
				if ($finalAlignment eq "") {
				  $finalAlignment .= ">$AncNum2AncName{$name}\n";
				}else{
				  $finalAlignment .= "\n>$AncNum2AncName{$name}\n";
				} 
			# otherwise, we have a species and not an ancestor
			}else{
				# always add a newline in front of a new fasta header, unless this is the very first 
				if ($finalAlignment eq "") {
				  $finalAlignment .= "$ali[$i]\n";
				}else{
				  $finalAlignment .= "\n$ali[$i]\n";
				} 
			}
		}else{
         $finalAlignment .= $ali[$i];
		}
	}

	# remove the first alignment column as this was artifically added to add the species with complete deletions
	if ($removeFirstAlignmentColumn == 1) {
		my $aliWithoutFirstColumn = "";
		foreach my $line (split(/\n/, $finalAlignment)) {
			if ($line =~ /^>/) {
				$aliWithoutFirstColumn .= "$line\n";
			}else{
				$aliWithoutFirstColumn .= substr($line, 1,length($line)) . "\n";
			}
		}
		print "\nFinal ali and anc WITH THE ARTIFICIAL alignment column:\n$finalAlignment\n" if ($verbose);
		$finalAlignment = $aliWithoutFirstColumn;
	}
	
	print "\nFinal ali and anc:\n$finalAlignment\n" if ($verbose);

	# write to BDB; 1 means use lockfile
	print "\nWrite result to $BDBFile with the key: $ElementID ... ";
	secureWriteBDB($BDBFile, $ElementID, $finalAlignment, 1, ".");
	print "DONE\n";

	# cleanup, unless $keepTemporaryFiles is set
	`rm -f ${tmpFile}*` if (! $keepTemporaryFiles);
	print "*** ALL DONE ***\n";

}

#######################################################################
# add the species (and their ancestors in case of sister species) that have a complete deletion
# For such species, we have no sequence to give to prank, therefore they are not in the output. 
# --> add those species back
# We do this by adding each complete-deletion species to the alignment and then run prank again to estimate the ancestors
# Because you cannot add an empty sequence, we add a new first alignment column that is an 'A' in all species with sequence or with a complete deletion
# This first column is removed after estimating the ancestors
#######################################################################
sub addCompleteDeletions {
	my $tmpFile = shift;

	my @completeDeletionList = keys %completeDeletionSpecies;
	$" = ","; 
	print "\nAdd these species with complete deletions: @completeDeletionList\n";
	# nothing to do, if no species has a complete deletion
	return 0 if ($#completeDeletionList == -1);

	# read the prank alignment of the species (no ancestors)
	# get the alignment length by just reading the aligned sequences
	my $call = "set -o pipefail; ConvertFastaToFastaOneLine.perl ${tmpFile}.best.fas | grep \">\" -v";
	print "$call\n" if ($verbose);
	my @ali = `$call`;
	cleanDie("ERROR: $call failed\n", $tmpFile) if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);	
	my $alignmentLength = -1;	
	chomp(@ali);
	for (my $i=0; $i<=$#ali; $i++) {
		cleanDie("ERROR: The alignment length of node differs from the length of other alignments: ", length($ali[$i]), " vs. $alignmentLength  ($ali[$i])\n", $tmpFile) if ( ($alignmentLength != -1) && (length($ali[$i]) != $alignmentLength) );
		$alignmentLength = length($ali[$i]);
	}
	print "\nalignmentLength: $alignmentLength\n" if ($verbose);

	my $completeDeletionAlignment = "A" . "-" x $alignmentLength;

	# now add the complete deletion species and the new alignment column to the alignment
	# use $tmpFile.addDelSpecies
	open FILEIN, "${tmpFile}.best.fas" || cleanDie("ERROR: cannot read from ${tmpFile}.best.fas\n", $tmpFile);
	@ali = <FILEIN>;
	close FILEIN;
	chomp(@ali);
	open FILEOUT, ">${tmpFile}.addDelSpecies.best.fas" || cleanDie("ERROR: cannot write to ${tmpFile}.addDelSpecies.best.fas\n", $tmpFile);
	my $currentName = "";		# proper name of the species or internal node
	for (my $i=0; $i<=$#ali; $i++) {
		# fasta header vs sequence
		if ($ali[$i] =~ /^>/) {
			print FILEOUT "$ali[$i]\nA";		# add the new alignment column
		}else{
			print FILEOUT "$ali[$i]\n";
		}
	}
	# now add the deletion species
	for my $species (@completeDeletionList) {
		print FILEOUT ">$species\n$completeDeletionAlignment\n";
	}
	close FILEOUT;
	
	`cp ${tmpFile}.addDelSpecies.best.fas ${tmpFile}.addDelSpecies.best.fasXXX`;
	
	# now replace the original alignment without the deletion species with the new alignment file with those species
	$call = "set -o pipefail; mv ${tmpFile}.addDelSpecies.best.fas ${tmpFile}.best.fas";
	print "$call\n" if ($verbose);
	system("$call") == 0 || cleanDie("ERROR: $call failed\n", $tmpFile);
	
	# this will call prank and only estimate ancestors (no realignment), output files will be $tmpFile.best.anc.fas etc
	# NOTE: There is a random element in prank of where to place ambiguous indels. Therefore, lets fix the seed so that we always get the same.
	$call = "set -o pipefail; prank -d=${tmpFile}.best.fas -t=${tmpFile}.tree.nh -o=${tmpFile}.best -keep -showtree -showanc -prunetree -seed=10";
	runPrankAndCheckForError($call, $tmpFile);

	return 1;
}
