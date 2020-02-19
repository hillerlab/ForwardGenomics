#!/usr/bin/perl

package MyBDB;
use strict;
use Exporter;
use Globals;
use BerkeleyDB;

our @ISA = ('Exporter');
our @EXPORT = qw(secureWriteBDB readBDB);


##################################################################
# get value from a DBD file for a given key
##################################################################
sub readBDB {
	my ($file, $key) = @_;

	my $Path = `pwd`;  chomp($Path);
	$key =~ s/$Path//g;		# get rid of the absolute path until chr*/CNE*/
	$key =~ s/\/\//\//g;		# get rid of // in e.g. chr2/CNE.2//CNE.2  
	$key = substr($key, 1, length($key)) if (substr($key, 0,1) eq "/");			# get rid of leading /

	my %h;
	tie %h, "BerkeleyDB::Btree", -Filename => "$file", -Flags => DB_RDONLY or die "Cannot open BerkeleyDB: $file   $! $BerkeleyDB::Error\n";
	my $value = "";
	if (exists $h{$key}) {
		$value = $h{$key};
	}
	untie %h;
	return $value;
}


##################################################################
# lock the DBD file, then write and unlock
##################################################################
sub secureWriteBDB {
	my ($file, $key, $value, $useLock, $dir) = @_;

	my $Path = `pwd`;  chomp($Path);
	$key =~ s/$Path//g;		# get rid of the absolute path until chr*/CNE*/
	$key =~ s/\/\//\//g;		# get rid of // in e.g. chr2/CNE.2//CNE.2  
	$key = substr($key, 1, length($key)) if (substr($key, 0,1) eq "/");			# get rid of leading /
	
#	my $filename = `basename $file`;
#	if ($useLock) {
#		print "lockfile -1 $dir/lockFile.$filename\n";
#		system "lockfile -1 $dir/lockFile.$filename" || die "ERROR: cannot get lock file: $dir/lockFile.$filename";
#	}
	my %hash;
	tie %hash, "BerkeleyDB::Btree",	-Filename => "$file", -Flags => DB_CREATE or die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";
	$hash{$key} = $value;
	untie %hash;

	# for lustre. Make sure it writes everything to the disk.
	# Ju found on April 30, 2017 that jobs stall if sync is called. Therefore, we skip that step now.
	#system("sync");

	# now read the same key from the BDB and see if we get exactly the same. Otherwise die. 
	my %h;
	tie %h, "BerkeleyDB::Btree", -Filename => "$file", -Flags => DB_RDONLY or die "Cannot open BerkeleyDB: $file   $! $BerkeleyDB::Error\n";
	my $valueRead = "";
	if (exists $h{$key}) {
		$valueRead = $h{$key};
	}
	untie %h;
	die "ERROR in secureWriteBDB: reading $key from $file does not give the value that was stored\nINPUT $value\nOUTPUT $valueRead\n" if ($value ne $valueRead);
        
	# Ju found on April 30, 2017 that jobs stall if sync is called. Therefore, we skip that step now.
        #system("sync");

#	if ($useLock) {
#		print "rm -f $dir/lockFile.$filename\n";
#		system "rm -f $dir/lockFile.$filename" || die "ERROR: cannot delete $dir/lockFile.$filename";
#	}
}


=pod
Should not be used anymore. Use secureWriteBDB above
##################################################################
# lock the DBD file, then write the entire hash and unlock
# $useLock is a flag: if 1 unix lockfile mechanism is used
##################################################################
sub secureWriteBDBHash {
	my ($file, $InputHash, $useLock, $dir) = @_;

	my $Path = `pwd`;  chomp($Path);
	
	if ($useLock) {
		system "lockfile -1 $dir/lockFile.$file" || die "ERROR: cannot get lock file: $dir/lockFile.$file";
	}
	my %hash;
	tie %hash, "BerkeleyDB::Btree",	-Filename => "$file", -Flags => DB_CREATE or die "Cannot open BerkeleyDB: $! $BerkeleyDB::Error\n";
	
	foreach my $key (keys %{$InputHash}) {	
		my $value = $InputHash->{$key};
		$key =~ s/$Path//g;		# get rid of the absolute path until chr*/CNE*/
		$key =~ s/\/\//\//g;		# get rid of // in e.g. chr2/CNE.2//CNE.2  
		$key = substr($key, 1, length($key)) if (substr($key, 0,1) eq "/");			# get rid of leading /
		$hash{$key} = $value;
	}
	untie %hash;

	if ($useLock) {
		system "rm -f $dir/lockFile.$file" || die "ERROR: cannot delete $dir/lockFile.$file";
	}
}
=cut

1;
