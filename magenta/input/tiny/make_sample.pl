#!/usr/bin/env perl
use 5.012;
use warnings;
use BioX::Seq::Stream;
use File::Basename;
use File::Spec;
my $read_length = 150;
my $frag_length = 500;

my ($recipe_file, $source_dir, $opt3) = @ARGV;

$source_dir = File::Spec->rel2abs( $source_dir);
die "Missing arguments: <recipe-file> <source-dir>\n" unless (-d "$source_dir");

my $out_base = $opt3 // $recipe_file . '_';	
open(my $I, '<', "$recipe_file") || die " Unable to read $recipe_file.\n";

while (my $line = readline($I) ) {
	chomp($line);
	next if ($line =~/^#/);
	my ($genome, $coverage, $from, $to) = split /\t/, $line;
	die "Wrong input format: expecting 4 columns in $recipe_file\n" unless ($to);
	my $chromosome = get_chr($genome);
	
	my $miniref =   $chromosome->range($from,$to)->seq;
	$chromosome->id = $genome;

	my $total_bases = length($miniref) * $coverage;
	my $printed_bases = 0;
	my $c = 0;
	say STDERR ">${genome}_${from}_${to}\tsize=", length($miniref), "\n$miniref";
	while ($printed_bases  <= $total_bases) {
		$c++;
		my $pos = int(rand(length($miniref) - $frag_length - 1));

		$printed_bases += 2 * $read_length;
		my $for = substr($miniref, $pos, $read_length);
		my $rev = rc(substr(substr($miniref, $pos, $frag_length), -1 * $read_length));
		say '@', "${genome}:$from:$to:$c/1", "\n",
			$for, "\n+\n!", 'I' x ($read_length - 1);
		say '@', "${genome}:$from:$to:$c/2", "\n",
			$rev, "\n+\n!", 'I' x ($read_length - 1);	

	}
}
sub get_chr{
	my ($genome) = @_;
	my @files = <$source_dir/$genome/*.fna>;
	die "[$genome] FASTA file not found at $source_dir/$genome/*.fna" unless $files[0];
	my $parser = BioX::Seq::Stream->new( $files[0] );
	return $parser->next_seq;
}
sub rc {

	my $seq = reverse(uc($_[0]));
	die unless $seq;
	$seq =~tr/ACGT/TGCA/;
	return $seq;
}