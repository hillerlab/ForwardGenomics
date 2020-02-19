#!/usr/bin/perl

# Michael Hiller, 2009

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

sub usage {
	die "usage: $0 inputFile\nall assembly names are replaced by species names\n";
};
my $NoSpaces = 0;
# options
GetOptions ("NoSpaces" => \$NoSpaces) || usage();
usage() if ($#ARGV < 0);

open(file1, $ARGV[0]);
my @Lines = <file1>;
close file1;

my $lines = "@Lines";

if ($NoSpaces == 0) {
$lines =~ s/anoCar1/lizard     /g;
$lines =~ s/bosTau3/cow        /g;
$lines =~ s/bosTau4/cow        /g;
$lines =~ s/calJac1/marmoset   /g;
$lines =~ s/canFam2/dog        /g;
$lines =~ s/cavPor2/guineaPig  /g;
$lines =~ s/cavPor3/guineaPig  /g;
$lines =~ s/choHof1/sloth      /g;
$lines =~ s/danRer4/zebrafish  /g;
$lines =~ s/danRer5/zebrafish  /g;
$lines =~ s/dasNov1/armadillo  /g;
$lines =~ s/dasNov2/armadillo  /g;
$lines =~ s/dipOrd1/kangarooRat/g;
$lines =~ s/echTel1/tenrec     /g;
$lines =~ s/equCab1/horse      /g;
$lines =~ s/equCab2/horse      /g;
$lines =~ s/eriEur1/hedgehog   /g;
$lines =~ s/felCat3/cat        /g;
$lines =~ s/galGal3/chicken    /g;
$lines =~ s/gasAcu1/stickleback/g;
$lines =~ s/gorGor1/gorilla    /g;
$lines =~ s/loxAfr1/elephant   /g;
$lines =~ s/loxAfr2/elephant   /g;
$lines =~ s/micMur1/mouseLemur /g;
$lines =~ s/monDom4/opossum    /g;
$lines =~ s/monDom5/opossum    /g;
$lines =~ s/myoLuc1/microbat   /g;
$lines =~ s/ochPri2/pika       /g;
$lines =~ s/ornAna1/platypus   /g;
$lines =~ s/oryCun1/rabbit     /g;
$lines =~ s/oryLat1/medaka     /g;
$lines =~ s/oryLat2/medaka     /g;
$lines =~ s/otoGar1/bushbaby   /g;
$lines =~ s/panTro2/chimpanzee /g;
$lines =~ s/petMar1/lamprey    /g;
$lines =~ s/ponAbe2/orangutan  /g;
$lines =~ s/proCap1/rockHyrax  /g;
$lines =~ s/pteVam1/megabat    /g;
$lines =~ s/rheMac2/rhesus     /g;
$lines =~ s/sorAra1/shrew      /g;
$lines =~ s/speTri1/squirrel   /g;
$lines =~ s/taeGut1/zebraFinch /g;
$lines =~ s/tarSyr1/tarsier    /g;
$lines =~ s/tetNig1/tetraodon  /g;
$lines =~ s/tupBel1/treeShrew  /g;
$lines =~ s/turTru1/dolphin    /g;
$lines =~ s/vicPac1/alpaca     /g;
$lines =~ s/xenTro2/frog       /g;
$lines =~     s/fr2/fugu       /g;
$lines =~     s/rn4/rat        /g;
$lines =~     s/rn5/rat        /g;
$lines =~     s/mm8/mouse      /g;
$lines =~     s/mm9/mouse      /g;
$lines =~    s/hg19/humanHg19  /g;
$lines =~    s/hg18/humanHg18  /g;

}else{

$lines =~ s/anoCar1/lizard/g;
$lines =~ s/bosTau3/cow/g;
$lines =~ s/bosTau4/cow/g;
$lines =~ s/calJac1/marmoset/g;
$lines =~ s/canFam2/dog/g;
$lines =~ s/cavPor2/guineaPig/g;
$lines =~ s/cavPor3/guineaPig/g;
$lines =~ s/choHof1/sloth/g;
$lines =~ s/danRer4/zebrafish/g;
$lines =~ s/danRer5/zebrafish/g;
$lines =~ s/dasNov1/armadillo/g;
$lines =~ s/dasNov2/armadillo/g;
$lines =~ s/dipOrd1/kangarooRat/g;
$lines =~ s/echTel1/tenrec/g;
$lines =~ s/equCab1/horse/g;
$lines =~ s/equCab2/horse/g;
$lines =~ s/eriEur1/hedgehog/g;
$lines =~ s/felCat3/cat/g;
$lines =~ s/galGal3/chicken/g;
$lines =~ s/gasAcu1/stickleback/g;
$lines =~ s/gorGor1/gorilla/g;
$lines =~ s/loxAfr1/elephant/g;
$lines =~ s/loxAfr2/elephant/g;
$lines =~ s/micMur1/mouseLemur/g;
$lines =~ s/monDom4/opossum/g;
$lines =~ s/monDom5/opossum/g;
$lines =~ s/myoLuc1/microbat/g;
$lines =~ s/ochPri2/pika/g;
$lines =~ s/ornAna1/platypus/g;
$lines =~ s/oryCun1/rabbit/g;
$lines =~ s/oryLat1/medaka/g;
$lines =~ s/oryLat2/medaka/g;
$lines =~ s/otoGar1/bushbaby/g;
$lines =~ s/panTro2/chimpanzee/g;
$lines =~ s/petMar1/lamprey/g;
$lines =~ s/ponAbe2/orangutan/g;
$lines =~ s/proCap1/rockHyrax/g;
$lines =~ s/pteVam1/megabat/g;
$lines =~ s/rheMac2/rhesus/g;
$lines =~ s/sorAra1/shrew/g;
$lines =~ s/speTri1/squirrel/g;
$lines =~ s/taeGut1/zebraFinch/g;
$lines =~ s/tarSyr1/tarsier/g;
$lines =~ s/tetNig1/tetraodon/g;
$lines =~ s/tupBel1/treeShrew/g;
$lines =~ s/turTru1/dolphin/g;
$lines =~ s/vicPac1/alpaca/g;
$lines =~ s/xenTro2/frog/g;
$lines =~     s/fr2/fugu/g;
$lines =~     s/rn4/rat/g;
$lines =~     s/rn5/rat/g;
$lines =~     s/mm8/mouse/g;
$lines =~     s/mm9/mouse/g;
$lines =~    s/hg19/humanHg19/g;
$lines =~    s/hg18/humanHg18/g;
}


print "$lines\n";
