#!/usr/bin/perl

package Globals;
use Exporter;
use Math::BigInt;
use Math::BigFloat;
use Math::Complex;
use Math::Trig;
our @ISA = ('Exporter');

#my $genomePath = $ENV{'genomePath'};

# globals
our $verbose = 0;
#our $CLUSTALPath = "$genomePath/bin/x86_64/";
#our $CLUSTALBinary = $CLUSTALPath . "/clustalw2";


%GeneticCode = qw (TTT Phe TTC Phe TTA Leu TTG Leu TCT Ser TCC Ser TCA Ser TCG Ser TAT Tyr TAC Tyr TAA Ter TAG Ter TGT Cys TGC Cys TGA Ter TGG Trp CTT Leu CTC Leu CTA Leu CTG Leu CCT Pro CCC Pro CCA Pro CCG Pro CAT His CAC His CAA Gln CAG Gln CGT Arg CGC Arg CGA Arg CGG Arg ATT Ile ATC Ile ATA Ile ATG Met ACT Thr ACC Thr ACA Thr ACG Thr AAT Asn AAC Asn AAA Lys AAG Lys AGT Ser AGC Ser AGA Arg AGG Arg GTT Val GTC Val GTA Val GTG Val GCT Ala GCC Ala GCA Ala GCG Ala GAT Asp GAC Asp GAA Glu GAG Glu GGT Gly GGC Gly GGA Gly GGG Gly);

our @EXPORT = qw($verbose %GeneticCode $CLUSTALPath $CLUSTALBinary);
