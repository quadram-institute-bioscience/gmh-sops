#!/usr/bin/env perl
use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

# Defaults
my $genus = 'Unknown';
my $species = 'sp.'; 
my $min_ratio  = 60.0;


my $VERSION    = 0.1.005;
my $reportFile = shift @ARGV;

my $top_species = `grep -w 'S' "$reportFile"  | head -n 1`;
$top_species =~s/^\s*//;

# Check top "S" species level
my ($ratio, undef, undef, undef, undef, @classification) = split /\s+/, $top_species;
if (defined $ratio and $ratio >= $min_ratio) {
    $genus = shift @classification;
    $species = shift @classification;
    # Genus as "E." rather than "Escherichia"
    $genus = uc( substr($genus, 0, 1) ) . '.';
} else {
    # Check top "G" genus level, eg:
    # 83.69  206578  7936    G       1578                    Lactobacillus
    my $top_species = `grep -w 'G' "$reportFile"  | head -n 1`;
    $top_species =~s/^\s*//;
    my ($ratio, undef, undef, undef, undef, @classification) = split /\s+/, $top_species;
    if (defined $ratio and $ratio >= $min_ratio) {
        $genus = shift @classification;
        $genus = uc( substr($genus, 0, 1) ) . '.';
    }
}



print " --genus $genus --species $species ";
