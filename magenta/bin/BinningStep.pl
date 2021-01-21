#!/usr/bin/env perl
use 5.012;
use warnings;
use Cwd qw(cwd);
use File::Basename;
use Data::Dumper;

my $job_dir = cwd;

my ($sampleID, $threads) = @ARGV;
die "Missing parameter: Sample_ID\n" unless ($sampleID);

$threads //= $ENV{'THREADS'};
$threads //= 1;

## INPUT
# Output of MetaBat (*.fa) and MaxBin (*.fasta)
## OUTPUT
# das_tool: 
#rm *unbinned.fa *tooShort.fa *lowDepth.fa
removeIfContains('unbinned.fa', 'tooShort.fa', 'lowDepth.fa');

### ==== RETRIEVE INPUT FILES FROM MAXBIN AND METABAT
# Metabat output (*.fa)
my @metabat_files = listByExt('fa');
# MaxBin output (*.fasta)
my @maxbin_files = listByExt('fasta');

say STDERR 'METABAT FILES: ', scalar @metabat_files, ":\n ", join(',', @metabat_files);
say STDERR 'MAXBIN FILES:  ', scalar @maxbin_files,  ":\n ", join(',', @maxbin_files);

### ==== PREPROCESS 
my $MaxBinCmd = qq(Fasta_to_Scaffolds2Bin.sh   --input_folder ./     --extension fasta > maxbin.tsv);
my $MetaBatCmd = qq(Fasta_to_Scaffolds2Bin.sh  --input_folder ./     --extension fa    > metabat.tsv);

if (scalar @maxbin_files) {
    say STDERR "INFO: Preprocessing MaxBin\n$MaxBinCmd";
    my $exit = system($MaxBinCmd);
    if (fileEmpty('maxbin.tsv')) {
        say STDERR "ERROR: MaxBin Preprocessing: no output or program execution failed: erasing array.";
        @maxbin_files=();
    }
}

if (scalar @metabat_files) {
    say STDERR "INFO: Preprocessing MaxBin:\n$MetaBatCmd";
    my $exit = system($MetaBatCmd);
    if (fileEmpty('metabat.tsv')) {
        say STDERR "ERROR: MaxBin Preprocessing: no output or program execution failed: erasing array.";
        @metabat_files=();
    } else {
        my $sedCmd = "sed -i 's/.metabat./_/g' metabat.tsv";
        system($sedCmd) || say STDERR "# Warning: 'sed' didn't work on metabat.tsv.";
    }
}
### ==== IF MAXBIN __AND__ METABAT:
my $labels = '';
my $tsv_files = '';

if (scalar @metabat_files and scalar @maxbin_files) {
    say STDERR "INFO: OK: Both MaxBin and MetaBat produced bins";
    $labels    = 'MaxBin2,MetaBat2';
    $tsv_files = 'maxbin.tsv,metabat.tsv';
} elsif (scalar @metabat_files) {
    say STDERR "INFO: Only MetaBat produced output";
    $labels    = 'MetaBat2';
    $tsv_files = 'metabat.tsv';
} elsif (scalar @maxbin_files) {
    say STDERR "INFO: Only MaxBin produced output";
    $labels    = 'MaxBin2';
    $tsv_files = 'maxbin.tsv';
} else {
    say STDERR "ERROR: Neither MaxBin nor MetaBat produced a valid output!";
    exit 1;
}

### ====EXECUTE DAS TOOL
my $das_command = qq(DAS_Tool -l $labels -i $tsv_files ) .
                  qq( -c contigs.fna  -o refine  --search_engine diamond  -t $threads  --write_bins);

system($das_command) || say STDERR "Warning: DAS_Tool returned non 0 (ignoring)";

### CHECK OUTPUT FILES

my @das_bins = listFiles("$job_dir/refine_DASTool_bins");

say STDERR scalar @das_bins, ' bins: ', join(',', @das_bins);

 
sub fileEmpty {
    my $filename = shift @_;
    if ( -e "$filename") {
        if ( -s "$filename") {
            return 0;
        } else {
            say STDERR "# Warning: file <$filename> not is empty.";
            return 1;
        }
    } else{
        say STDERR "# Warning: file <$filename> not found.";        
    }
    
}
sub listByExt {
    my ($ext, $directory) = @_;
    $directory //= $job_dir;
    my @files = listFiles($directory);
    my @filtered = ();
    for my $file (@files) {
        push(@filtered, $file) if ($file =~/\.$ext$/);
    }
    return @filtered;
}

sub removeIfContains { 
    my @files = listFiles();
    for my $file (@files) {
        for my $string (@_) {
            unlink("$job_dir/$file") if ($file =~/$string/);
        }
    }
}

sub listFiles {
    my $input_dir = $_[0] // $job_dir;
    opendir my $dir, "$job_dir" or die "ERROR: Cannot open directory <$job_dir>: $!";
    my @files = readdir $dir;
    closedir $dir;
    return @files;
}

sub checkBin {
    my ($cmd, $expected_status) = @_;
    $expected_status // 0;
    my $output = `$cmd`;
    my $exitstatus = $?;

}

__END__
  Fasta_to_Scaffolds2Bin.sh  \
          --input_folder ./  \
          --extension fa > metabat.tsv

  Fasta_to_Scaffolds2Bin.sh  \
          --input_folder ./  \
          --extension fasta > maxbin.tsv

  sed -i 's/.metabat./_/g' metabat.tsv
  DAS_Tool \
          -l MaxBin2,MetaBat2 \
          -i maxbin.tsv,metabat.tsv \
          -c contigs.fna \
          -o refine \
          --search_engine diamond \
          -t ${task.cpus} \
          --write_bins

  # OUTPUT MANAGEMENT
  # Dir:  refine_DASTool_bins/*.fa
  COUNTER=0
  for BIN in refine_DASTool_bins/*.fa*;
  do
    COUNTER=\$((COUNTER+1))
    mv \$BIN ${name}.das.bin_\$COUNTER.fa

  done
  # PDFs: refine_DASTool_scores.pdf, refine_DASTool_hqBins.pdf
  # QC:   refine_DASTool_summary.txt
  mv refine_DASTool_scores.pdf  ${name}.das.scores.pdf
  mv refine_DASTool_hqBins.pdf  ${name}.das.qc.pdf
  mv refine_DASTool_summary.txt ${name}.das.summary.txt

  output:
  set val(name), file("*.summary.txt"), file("*.pdf") into dasqc_ch
  set val(name), file("*.das.bin*.fa") into das_ch
---  

DAS_Tool -i methodA.scaffolds2bin,...,methodN.scaffolds2bin
         -l methodA,...,methodN -c contigs.fa -o myOutput

   -i, --bins                 Comma separated list of tab separated scaffolds to bin tables.
   -c, --contigs              Contigs in fasta format.
   -o, --outputbasename       Basename of output files.
   -l, --labels               Comma separated list of binning prediction names. (optional)
   --search_engine            Engine used for single copy gene identification [blast/diamond/usearch].
                              (default: usearch)
   --write_bin_evals          Write evaluation for each input bin set [0/1]. (default: 1)
   --create_plots             Create binning performance plots [0/1]. (default: 1)
   --write_bins               Export bins as fasta files  [0/1]. (default: 0)
   --proteins                 Predicted proteins in prodigal fasta format (>scaffoldID_geneNo).
                              Gene prediction step will be skipped if given. (optional)
   --score_threshold          Score threshold until selection algorithm will keep selecting bins [0..1].
                              (default: 0.5)
   --duplicate_penalty        Penalty for duplicate single copy genes per bin (weight b).
                              Only change if you know what you're doing. [0..3]
                              (default: 0.6)
   --megabin_penalty          Penalty for megabins (weight c). Only change if you know what you're doing. [0..3]
                              (default: 0.5)
   --db_directory             Directory of single copy gene database. (default: install_dir/db)
   --resume                   Use existing predicted single copy gene files from a previous run [0/1]. (default: 0)
   --debug                    Write debug information to log file.
   -t, --threads              Number of threads to use. (default: 1)
   -v, --version              Print version number and exit.
   -h, --help                 Show this message.


  sub fileEmpty "
  Can't use global @_ in "my" at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 86, near "shift @_"
  Global symbol "$filename" requires explicit package name at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 87.
  Global symbol "$filename" requires explicit package name at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 88.
  Global symbol "$filename" requires explicit package name at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 91.
  Global symbol "$filename" requires explicit package name at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 95.
  syntax error at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 98, near "}"
  Can't use global @_ in "my" at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 100, near "= @_"
  Global symbol "$directory" requires explicit package name at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 101.
  syntax error at /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl line 108, near "}"
  /qib/platforms/Informatics/GMH/nextflow/magenta_dev/bin/BinningStep.pl has too many errors.
