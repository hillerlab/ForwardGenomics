#!/usr/bin/perl

# Michael Hiller, 2015
# revamp of the %id calculation.
# We now expect an alignment together with all reconstructed internal nodes and compute the global and local %ids.

use strict;
use warnings;
use lib "./";
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use MyBDB;
use Bio::SeqIO;
use Bio::Tree::Draw::Cladogram;
use Bio::Tree::Tree;
use Bio::Tree::NodeI;
use Bio::Tree::TreeI;
use Bio::Tree::TreeFunctionsI;
use Bio::TreeIO;
use POSIX;

# options
my $verbose = 0;
my $verbose1 = 0;

my $global = 0;		# flags
my $local = 0;

# low qual parameters
my $doNotIgnoreN = 0;				# if set to 1 do not ignore N-containing alignment columns but treat as a always-mismatch nt
my $maxNumLowQual = 50;				# ignore species if it has more than 50 low qual bases
my $maxFractionLowQual = 0.20; 	# ignore species if it has more than 20% low qual (counted only for pos where species and ancestor has a BASE)

# whether to cap insertions above a certain length
my $capInsLongerThan = INT_MAX;

# input, temp and output files
my $keepTemporaryFiles = 0;
my $GlobalBDBFile = "";
my $LocalBDBFile = "";
my $treeFile = "";
my $elementsFile = "";				

# print options
my $brief = 0;
my $displayAli = 0;					# provide a view of the sorted alis with the ancestor (. = identity)
my $LowQualOutputChar = "?";		# fill in the -display ali low quality bases with that char

# ancestors
my $allowedAncestralNodes = "";	# specifies which internal nodes can be used as the ingroup ancestor
my $requireNoOutgroup = 0;			# if set to 1, no outgroup is required and thus the root can be ingroup-ancestor

my $speciesWithoutCNE = "";		# if given that overwrites the speciesListSorted_GetPercentID file


sub usage {
	die "usage: $0 Fasta_or_BDBFile ElementID -treeFile string -allowedAncestralNodes string  [-global | -local] [PARAMETERS]
\
$0 will read the Prank-alignment from either a BDB file or a fasta file (must have the suffix .fa) produced by Maf2SpanningSeq_PRANK.perl
It computes global and local \%id values.
If a sequence is completely deleted in an ancestor, all downstream branches get \%id = 0. 	
\
general parameters:
	-treeFile string                          file with the tree in phastCons model or in newick format (e.g. tree.mod or tree.nh). Note: all ancestors must be full named (tree_doctor -a).
	-doNotIgnoreN                             do not ignore alignment columns with N (low quality) but count them as mismatches/indels
	-keepTemporaryFiles                       flag: keep the files (pruned tree)
	-elementsFile                             optional: bed file listing the elements coordinates. If given, we will output the length of the element as the 4th column for the local \%id values.
\
parameters for \%id reconstruction
	-local                                    compute local \%id values. Output goes to stdout unless you give -LocalBDBFile string
	-global                                   compute global \%id values. Output goes to stdout unless you give -GlobalBDBFile string
	-allowedAncestralNodes string             comma-separated list of node names that can be used as the ingroup ancestor
	-requireNoOutgroup                        flag: if set the root can also be the ingroup ancestor (in this case you don't need outgroup species)
	-capInsLongerThan int                     above this threshold, inserted bases will ignored (e.g. if set to 10, a 15 bp insertion will be counted as a 10 bp insertion) (default: $capInsLongerThan)
\
parameters for output
	-GlobalBDBFile string                     name of the BDBFile to output the global \%id values. If not given, output to stdout.
	-LocalBDBFile string                      name of the BDBFile to output the local \%id values. Must be different from GlobalBDBFile. If not given, output to stdout. 
\
print options
	-brief                                    no output if you set the -GlobalBDBFile and -LocalBDBFile
	-displayAli                               show a visual representation of the local and global ancestor-referenced alignment
	-v|verbose                                verbose output
	-v1|verbose1                              verbose the \%id calculation per alignment column
\
\n"};
GetOptions ("v|verbose"  => \$verbose, "v1|verbose1"  => \$verbose1, "keepTemporaryFiles" => \$keepTemporaryFiles, "treeFile=s" => \$treeFile,
				"GlobalBDBFile=s" => \$GlobalBDBFile, "LocalBDBFile=s" => \$LocalBDBFile, 
				"doNotIgnoreN" =>\$doNotIgnoreN, "displayAli" => \$displayAli, "elementsFile=s" => \$elementsFile,
			   "allowedAncestralNodes=s" => \$allowedAncestralNodes, "requireNoOutgroup" => \$requireNoOutgroup,
				"LowQualOutputChar=s" => \$LowQualOutputChar, "global"  => \$global, "local"  => \$local,
				"speciesWithoutCNE=s" => \$speciesWithoutCNE, "brief" => \$brief,
				"maxFractionLowQual=f" => \$maxFractionLowQual, "capInsLongerThan=i" => \$capInsLongerThan
				) || usage();
if ($#ARGV < 1) {
	print STDERR "You must give inputFile and ElementID\n\n"; 
	usage();
}
if ($treeFile eq ""){ 
	print STDERR "You must set -treeFile filename\n\n"; 
	usage();
}
if ($allowedAncestralNodes eq "") {
	print STDERR "You must set -allowedAncestralNodes string\n\n"; 
	usage();
}
die "ERROR: You must set either -global or -local\n" if ($global == 0 && $local == 0);
die "ERROR: you cannot set the global and local BDB file to the same name\n" if ($GlobalBDBFile eq $LocalBDBFile && $LocalBDBFile ne "");

if ($verbose) {
	$| = 1;		# == fflush(stdout)
}
my $ElementID = $ARGV[1];
my $curSpecies = "";		# global var just for verbose output purpose

# $ancestorHash is a pointer to a hash that contains all possible allowed ingroup-ancestral nodes
my @ancestorList = split(/,/, $allowedAncestralNodes);
my $ancestorHash = arrayToHash(@ancestorList);

# list of all leaves in the tree. Needed to have consistent output of all these species for the global %id format (now a table elements x species).
my @allLeaves; 
# this function fills the @allLeaves array
getListOfLeavesFromTree();

my $elementLength = getElementLength($elementsFile);

# now compute %id
doit($ARGV[0]);

#########################################################################################################
## given an array, the function inits a hash with each element (hash value set to 1)
## returns pointer to hash
##########################################################################################################
sub arrayToHash {
	my %H;
	foreach my $element (@_) {
		$H{$element} = 1;
	};
	return \%H;
}


##################################################################
# if the argument is a bed file name open it and extract the length of the element
# otherwise return -1
##################################################################
sub getElementLength {
	my $elementsFile = shift;
	return -1 if ($elementsFile eq "");

	open(file1, $elementsFile) || die "ERROR: cannot open $elementsFile\n";
	my $line1 = "";
	while ($line1 = <file1>) {
		chomp ($line1);
		my @a = split(/\t/, $line1);
		if ($a[3] eq $ElementID) {
			close(file1);
			return $a[2]-$a[1];
		}
	}
	close(file1);
	return -1;
}

##################################################################
# decide whether a given node name is a leave based on the presence of the char "-"
# mm10, rn5 are leaves,  mm10-rn5 is not a leave
##################################################################
sub isLeaveByString {
	my $node = shift;
	if (index($node, "-") != -1) {
		return 0;
	}else{
		return 1;
	}
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
# given the tree.mod or tree.nh, get a full list of all leaves
#######################################################################
sub getListOfLeavesFromTree {
	my $tmpFile = `mktemp /dev/shm/tree.nh.XXXXXXXXX`;
	chomp($tmpFile);

	my $call = "set -o pipefail; ConvertTreeModToNewick.perl $treeFile $tmpFile";
	print "$call\n" if ($verbose);
	system("$call") == 0 || cleanDie("ERROR: $call failed\n", $tmpFile);	

	# read the tree
	my $BioTree = Bio::TreeIO->new('-format' => 'newick', '-file' => $tmpFile);
	my $tree = $BioTree->next_tree;

	my @taxa = $tree->get_leaf_nodes;
	my @allLeavesTmp;
	for my $node ( @taxa ) {
   	push @allLeavesTmp, $node->id;
		print "All leaves in full tree: ", $node->id, "\n" if ($verbose);
	}
	
	# clean tmpFile
	`rm -f $tmpFile`;
	
	# sort the array --> allLeaves is a global var
	@allLeaves = sort{ $a cmp $b } @allLeavesTmp;
}		


##################################################################
# to a top-down tree traversal and return the first node that occurs in the $ancestorHash --> this is the oldest ancestor
# If this ancestor is also the root (we have no outgroup) only return this node if -requireNoOutgroup is set
##################################################################
sub findIngroupAncestralSeq {
	my $tree = shift;
	my $rootName = $tree->get_root_node->id;

	print "\nrun findIngroupAncestralSeq() given the following possible ingroup ancestors: $allowedAncestralNodes\n" if ($verbose);
	for my $node ( $tree->get_nodes(-order => 'depth')) {
		print "Visit node: ", $node->id, "\n" if ($verbose);
		if (exists ($ancestorHash->{$node->id})) {
			print "==> found oldest ingroup ancestor ", $node->id if ($verbose);
			# if this ancestor == root, only use it if -requireNoOutgroup is set
			if ($node->id eq $rootName) {
				print " However that is the root !! (-requireNoOutgroup $requireNoOutgroup)\n" if ($verbose);
				if ($requireNoOutgroup == 1) {
					return $node;
				}
			}else{
				print " (different than the root)\n" if ($verbose);
				return $node;
			}
		}
	}
	# if we cannot find the right ingroup-ancestor, we print an error message but exit(0). 
	print STDERR "ERROR: could not find a ingroup ancestor ", ($requireNoOutgroup == 1 ? "" : "that is different from the root (an proper outgroup is present) "), "for $ElementID\n";
	exit(0);
}

##################################################################
# given a tree and an ancestral node, fill a hash with descendants
# if $onlyLeaves == 1, consider only leave nodes
##################################################################
sub getIngroupDescendants {
	my ($tree, $ingroupAncestorNode, $onlyLeaves) = @_;
	my %H;

	# all descendents of this node
	print "\nGet each descendent ", ($onlyLeaves == 1 ? "(only leaves) " : ""), "of node: ", $ingroupAncestorNode->id, "\n" if ($verbose);
	for my $child ( $ingroupAncestorNode->get_all_Descendents ) {
		if ($onlyLeaves == 1 && ! $child->is_Leaf) {
			print "\t IGNORE non-leaf ", $child->id, "\n" if ($verbose);
			next;
		}
		print "\t", $child->id, "\n" if ($verbose);
		$H{$child->id} = 1;
	}
	return \%H;
}

##################################################################
# determine %id for all species for that file
##################################################################
sub doit {
	my $file = shift;

	# read the alignment from a file or BDB
	my @ali;
	if ($file =~ /\.fa$/) {
		open(file1, $file) || die "ERROR: cannot open $file\n";
		@ali = <file1>;
		close file1;
		chomp(@ali);
	}else{
		my $aliFromBDB = readBDB($file, $ElementID); 
		@ali = split(/\n/, $aliFromBDB);
	}

	if ($verbose) {
		$" = "\n";
		print "ali read from $file:\n@ali\n\n";
	}

	# save exit and no output at all for those elements that could not be aligned
	if ($ali[0] eq "NOTHING_TO_ALIGN") {
		print "no alignment in $ARGV[0] for $ElementID: $ali[0]\n" if (! $brief);
		exit(0);
	}

	# now read the alignment into a hash, double check if the input format is correct and all alignments are of the same length
	# also, get a list of species (leaves)
	my %node2Ali;
	my %leaves;			# hash of leaves
	my $alignmentLength = -1;	
	for (my $i=0; $i<=$#ali; $i+=2){
		my $header = $ali[$i];
		my $aliLine = $ali[$i+1];
		
		if ($header !~ />(.*)/) {
			die "ERROR: expect a valid fasta header at line $i but got $header\nNOTE: the input must be one-line fasta\n";
		}else{
			my $node = $1;
			die "ERROR: The alignment length of $node differs from the length of other alignments: ", length($aliLine), " vs. $alignmentLength\n" if ( ($alignmentLength != -1) && (length($aliLine) != $alignmentLength) );
			die "ERROR: The alignment length of $node is 0\n" if ( length($aliLine) == 0 );
			$alignmentLength = length($aliLine);
			$node2Ali{$node} = $aliLine;
			print "\t$node\t$aliLine\n" if ($verbose);
			if (isLeaveByString($node)) {
				$leaves{$node} = 1;
				print "\t$node is a leave\n" if ($verbose);
			}
		}
	}

	# now prune the tree using tree_doctor. Only get those leaves that are in our alignment. 
	$" = ",";
	my @listOfLeaves = keys %leaves;
	my $listOfLeavesString = "@listOfLeaves";
	print "==> list of leaves: $listOfLeavesString\n" if ($verbose);
	my $tmpFile = `mktemp /dev/shm/prunedTree.XXXXXXXXX`;
	chomp($tmpFile);
	my $call = "set -o pipefail; ConvertTreeModToNewick.perl $treeFile /dev/stdout | tree_doctor /dev/stdin -n -P $listOfLeavesString -a > $tmpFile";
	print "$call\n" if ($verbose);
	system("$call") == 0 || cleanDie("ERROR: $call failed\n", $tmpFile);	

	# read the pruned tree
	my $BioTree = Bio::TreeIO->new('-format' => 'newick', '-file' => $tmpFile);
	my $tree = $BioTree->next_tree;
	# clean tmpFile
	`rm -f $tmpFile`;
	
	# get the root (needed for local %id)
	my $rootName = $tree->get_root_node->id;
	printf "\nread tree: root $rootName\n" if (! $brief);
	# sanity check
	die "ERROR: do not find $rootName in the alignment\n" if (! exists $node2Ali{$rootName});

	# find the correct ingroup-ancestral seq
	my $ingroupAncestorNode = findIngroupAncestralSeq($tree);
	my $ingroupAncestorName = $ingroupAncestorNode->id;
	my $ancestralSeq = "";
	$ancestralSeq = $node2Ali{$ingroupAncestorName} if (exists $node2Ali{$ingroupAncestorName});
	die "ERROR: cannot get ancestral sequence of ancestor $ingroupAncestorName\n" if ($ancestralSeq eq "");


	#######################
	# global %ids
	if ($global == 1) {
		my %aliLinesGlobal;
		my %percentIDGlobal;
		printf "compute global %%id values using ancestor: $ingroupAncestorName\n" if (! $brief);

		# get a pointer to a hash listing all descendant species from $ingroupAncestorNode (only leaves --> last parameter == 1) 
		my $ingroupSpecies = getIngroupDescendants($tree, $ingroupAncestorNode, 1);

		foreach my $leave (keys %leaves) {
			if (! exists $ingroupSpecies->{$leave}) {
				print "\nignore species $leave because it is not an ingroup wrt $ingroupAncestorName\n" if ($verbose);
				next;
			}
			die "ERROR: do not find leave $leave in the alignment\n" if (! exists $node2Ali{$leave});
			my $seq = $node2Ali{$leave};
			die "ERROR: The alignment of $leave differs in length from $ingroupAncestorName\n$ancestralSeq\n$seq" if (length ($seq) != length($ancestralSeq));
			print "\ncompute global %%id for $leave ........ \n" if ($verbose);
			$curSpecies = "$leave";
			# now compute the %id value
			my ($M, $len, $percentID, $aliLine) = getPercentID($seq, $ancestralSeq);
			# store the value if it is not -1 (meaning too much low qual)
			# Note: in outputGlobalPerIDs() we output 'NA' for all species without a %id value. 
			$percentIDGlobal{$curSpecies} = $percentID if ($percentID != -1);
			$aliLinesGlobal{$curSpecies} = $aliLine;
		}
		# output, either to BDB or to stdout
		print "\noutput global \%id values ....... \n" if ($verbose);
		outputGlobalPerIDs(\%percentIDGlobal);

		# display the global anc-referenced alignment
		if ($displayAli) {
			displayAli(\%aliLinesGlobal, $ancestralSeq, $ingroupAncestorName, 0);
		}

		print "\n#####################################################################################################################################################################################\n" if ($verbose);
	}

	#######################
	# local %ids
	if ($local == 1) {
		my %aliLinesLocal;
		my %percentIDLocal;
		printf "compute local %%id values using ancestor: $ingroupAncestorName\n" if (! $brief);

		# get a pointer to a hash listing all descendant nodes from $ingroupAncestorNode (leaves + internal nodes --> last parameter == 0) 
		# Note: the $ingroupAncestorNode is not in this hash. Therefore in the for loop below, we take a node and its ancestor and compute the local %id
		my $ingroupNodes = getIngroupDescendants($tree, $ingroupAncestorNode, 0);

		# traverse the tree and compute all local %id values
		for my $node ( $tree->get_nodes(-order => 'depth')) {
			my $nodeName = $node->id;
			next if ($nodeName eq $rootName);
			if (! exists $ingroupNodes->{$nodeName}) {
				print "\nignore node $nodeName because it is not an ingroup wrt $ingroupAncestorName\n" if ($verbose);
				next;
			}
			my $parentNodeName = $node->ancestor->id;
			die "ERROR: do not find node $nodeName in the alignment\n" if (! exists $node2Ali{$nodeName});
			die "ERROR: do not find node $parentNodeName in the alignment\n" if (! exists $node2Ali{$parentNodeName});
			my $seq = $node2Ali{$nodeName};
			my $parentSeq = $node2Ali{$parentNodeName};
			die "ERROR: The alignment of $nodeName differs in length from $parentNodeName\n$parentSeq\n$seq" if (length ($seq) != length($parentSeq));
			print "\ncompute local %%id for branch $parentNodeName - $nodeName ........ \n" if ($verbose);
			$curSpecies = "$parentNodeName - $nodeName";
			my ($M, $len, $percentID, $aliLine) = getPercentID($seq, $parentSeq);
			# percentID == -1 means too much low qual --> we output these values as 'NA' 
			# Reason: For Forward Genomics, we need to know which leaves had genomic (high or low quality) sequence to prune the full tree like the Prank step did. 
			# NOTE: We use here the endpoint of the branch (nodeName) as the key
			$percentIDLocal{$nodeName} = $percentID;
			$aliLinesLocal{$curSpecies} = $aliLine;
		}
		print "\noutput local \%id values ....... \n" if ($verbose);
		outputLocalPerIDs(\%percentIDLocal);

		# display the local anc-referenced alignment
		if ($displayAli) {
			displayAli(\%aliLinesLocal, $ancestralSeq, $ingroupAncestorName, 1);
		}
		print "\n#####################################################################################################################################################################################\n" if ($verbose);
	}

}


####################################################################################################
# output the global perIDs
# format: elementID hg18 panTro4 ...
# --> all values in one line
# NOTE: We always output the header line too. This allows to later either analyze this element separately or to combine the values of many elements.
####################################################################################################
sub outputGlobalPerIDs {
	my $percentIDGlobal = shift;		# pointer to hash
	
	my $NA = "NA";		# no value (either no alignment given or too much low qual)
	my $header = "species";
	my $values = "$ElementID";
	# Note: this array must be sorted to assure a consistent species order
	for (my $i=0; $i<=$#allLeaves; $i++) {
		$header .= " $allLeaves[$i]";
		my $value = $NA;
		if (exists $percentIDGlobal->{$allLeaves[$i]}) {
			$value = $percentIDGlobal->{$allLeaves[$i]};
			print "--> $allLeaves[$i] has a \%id value: $value\n" if ($verbose);
		}else{
			print "--> $allLeaves[$i] has NO \%id value: $value\n" if ($verbose);
		}
		$values .= " $value";
	}

	my $result = "$header\n$values\n";
	if ($GlobalBDBFile ne "" ){
		print "Write \n\t$header\n\t$values\n to $GlobalBDBFile with the key: $ElementID ... " if ($verbose);
		secureWriteBDB($GlobalBDBFile, $ElementID, $result, 1, ".");
		print "DONE\n" if ($verbose);
	}else{
		print "Global \%id values for $ElementID:\n$result\n";
	}
}

####################################################################################################
# output the local perIDs
# format: 
#		branch id pid 
#		Alpaca chr14.100 0.986842105263158
#		Alpaca-Dolphin chr14.100 1
# NOTE: We output the header line only if the output is to stdout
####################################################################################################
sub outputLocalPerIDs {
	my $percentIDLocal = shift;
	
	my $header = "branch id pid len";
	my $values = "";
	# Note: this array must be sorted to assure a consistent species order
	foreach my $branch (keys %{$percentIDLocal}) {
		if ($percentIDLocal->{$branch} == -1){ 
			$values .= "$branch $ElementID NA $elementLength\n";
		}else{
			$values .= "$branch $ElementID $percentIDLocal->{$branch} $elementLength\n"
		}
	}

	my $result = "$header\n$values";
	if ($LocalBDBFile ne "" ){
		print "Write \n\t$header\n\t$values\n to $LocalBDBFile with the key: $ElementID ... " if ($verbose);
		secureWriteBDB($LocalBDBFile, $ElementID, $result, 1, ".");
		print "DONE\n" if ($verbose);
	}else{
		print "Local \%id values for $ElementID:\n$result\n";
	}
}


####################################################################################################
# compare a seq to the ancestor seq (might contain gaps in both seqs at one pos) and return matches, ali length and percentID
# if there are too many low quality chars, we set percentID to -1
####################################################################################################
sub getPercentID {
	my ($seq, $ancestralSeq) = @_;

	my $M = 0; 								# matches
	my $len = 0;							# length of the alignment excluding positions that are both gaps
	my $numN = 0;							# number of low quality chars in the sequence
	my $numGapsInBoth = 0;				# count how many columns are - in both ancestor and sequence --> we need this for a special case in the local %id
	my $numRealBasesWoGaps = 0;		# number of real bases in the sequence that do not align to a - in the ancestor (match or mismatch)
	my $currentInsertionLength = 0;	# keeps track of the length of the current insertion. If longer than $capInsLongerThan, do not add the additional ins bases to $len
	my $aliLine = "";

	print "getPercentID\nseq:\t$seq\nanc:\t$ancestralSeq\n" if ($verbose);

	print "\tPOS   seq vs anc\n" if ($verbose1);
	for (my $i=0; $i<length($seq); $i++) {
		my $s = uc substr($seq, $i,1);
		my $a = uc substr($ancestralSeq, $i,1);
		print "\tPOS $i  $s vs $a  " if ($verbose1);

		# ignore column: no base in both the anc and seq
		if ($a eq "-" && $s eq "-") {		
			$aliLine .= " ";
			print "IGNORE double gap\n" if ($verbose1);
			$numGapsInBoth++;
			next;
		}

		# low quality base in the sequence
		$numN ++ if ($s eq "N");
		if ($s eq "N" && $doNotIgnoreN == 0) {
			if ($a ne "-") {		
				$aliLine .= $LowQualOutputChar;		# write that char if there is an ancestral base
			}else{
				$aliLine .= " ";
			}
			print "IGNORE low qual\n" if ($verbose1);
			next; 
		}

		# count the number of ancestral bases
#		$numAncestralBases ++ if ($a ne "-");

		# reset variable if the ancestor has a base
		$currentInsertionLength = 0 if ($a ne "-");

		# alignment length (NOTE: we skip over columns with - in seq and anc already above)
		$len ++;

		# count at how many positions there is a real base in the sequence and in the ancestor 
		# Note: We ignore here the columns with gaps, because in low qual regions insertions/deletions are more likely and some "real bases" could just be a bit above the quality threshold
		$numRealBasesWoGaps ++ if (($s =~ /[ACGT]/ && $a ne "-"));
		
		# Match or Mismatch/Indel
		# Match
		if ($a eq $s) {
			$M++;
			$aliLine .= ".";
			print "Match\n" if ($verbose1);
		# Insertion
		}elsif ($a eq "-" && $s ne "-") {
			$currentInsertionLength ++;
			if ($currentInsertionLength > $capInsLongerThan) {
				print "Ins longer than $capInsLongerThan ($currentInsertionLength bp inserted) \n" if ($verbose1);
				# We don't count this insertion --> do not add this alignment column to $len (because we already have, subtract)
				$len--;				
				$aliLine .= lc "$s";
			}else{
				print "Ins\n" if ($verbose1);
				$aliLine .= "$s";
			}	
		# Deletion or MisMatch
		}else {
			$aliLine .= "$s";
			print "Del or Mismatch\n" if ($verbose1);
		}		
	}
	my $fractionNChars = $numN/(   ($numRealBasesWoGaps+$numN) > 0 ? ($numRealBasesWoGaps+$numN) : 1   );
	my $percentID = -1;
	$percentID = $M/$len if ($len > 0);

	# special case for local %id:
	#  If a complete deletion happens in an internal node, then all downstream branches have just - for the ancestor and species. 
	#  In this case, set %id to 0 meaning "still a complete deletion"
	if (length($seq) == $numGapsInBoth) {
		$percentID = 0;
		print "SPECIAL CASE for local \%id: internal deletion: the sequence of this branch was already deleted in the ancestor --> set \%id = 0\n" if ($verbose);
		$aliLine = "-" x $numGapsInBoth;
	}

	# verbose output
	print "len:             $len\n" if ($verbose);
	print "matches:         $M\n" if ($verbose);	
	print "#N's:            $numN  maxNumLowQual: $maxNumLowQual  numN/(numRealBasesWoGaps+numN) = $numN/($numRealBasesWoGaps+$numN) = $fractionNChars\n" if ($verbose);
	print "percentID:       $percentID\n" if ($verbose);

	# check low quality:  no %id value if there is too much low quality 
	if ($numN > $maxNumLowQual || $fractionNChars > $maxFractionLowQual) {
		print "$ElementID $curSpecies  LowQual Flag: numNChars: $numN (max $maxNumLowQual) || fractionNChars: $fractionNChars (max $maxFractionLowQual)  --> set \%id to -1\n" if ($verbose);
		$M = 0; $len = 0; $percentID = -1;
		return ($M, $len, $percentID, $aliLine);
	}
	return ($M, $len, $percentID, $aliLine);

}



####################################################################################################
# displayAli
####################################################################################################
sub displayAli {
	my ($aliLines, $ancestralSeq, $ancestorName, $local) = @_; 

	print "alignment to the ancestral reconstructed sequence: $ElementID\n";
	print "LEGEND: \n\tThe first line gives the reconstructed ancestral sequence. \n";
	print "\tThen the sequences for the species are shown:\n";
	print "\t\t . means identical nucleotide wrt to ancestral sequence\n";
	print "\t\t [ACGT-] (upper case) means substitution or insertion/deletion wrt to ancestral sequence\n";
	print "\t\t [acgt]  (lower case) is are inserted bases above capInsLongerThan there were ignored\n";
	print "\t\t ? means low quality nucleotide. These positions are ignored.\n";
	
	my $spaces = 19;
	$spaces = 41 if ($local == 1);

	printf("%${spaces}s $ancestralSeq\n", $ancestorName);
	foreach my $key (keys %${aliLines}){
		printf("%${spaces}s $aliLines->{$key}\n", $key);
	}

=pod
	for (my $i=0; $i<=$#sorted; $i++) {
		if (exists $aliLines->{$sorted[$i]}) {
			my @s = `set -o pipefail; echo "$sorted[$i]" | Assembly2Species.perl /dev/stdin`;
			die "ERROR: echo $sorted[$i] | Assembly2Species.perl  FAILED\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
			chomp(@s);
			printf "%-16s %s\n",$s[0],$aliLines->{$sorted[$i]};
		}
	}
=cut
}

