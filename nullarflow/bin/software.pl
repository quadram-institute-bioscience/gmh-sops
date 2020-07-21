#!/usr/bin/env perl
# ABSTRACT - A tool to collect software versions, check software presence, print HTML/TXT report
use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

my $opt_html = 'versions_mqc.html';
GetOptions(
  'html=s'    => \$opt_html,
);
open(my $HTML, '>', "$opt_html") || die "Unable to write to $opt_html\n";

my $databases = {
    'pubmlst' => {
        desc  => 'PubMLST website (https://pubmlst.org/) sited at the University of Oxford.',
        uri   => 'https://pubmlst.org/',
        paper => ' Jolley & Maiden 2010, BMC Bioinformatics, 11:595',
    },
    'resfinder' => {
        desc => 'ResFinder database',
        uri  => 'doi:10.1093/jac/dks261',

    },
    'VFDB' => {
        desc => 'VFDB database',
        uri  => 'doi:10.1093/nar/gkv1239',
    }
};

my $software = {
    'roary' => {
        desc     => 'pangenome detection',
        cmd      => 'roary --version 2>/dev/null',
        regex    => '(\d+.\d+\.?\d+?)',
        paper    => 'Page AJ (2015)',
        doi      => 'http://dx.doi.org/10.1093/bioinformatics/btv421',
        uri      => 'https://sanger-pathogens.github.io/Roary/',
    },
    'fastp' => {
        cmd   => 'fastp --version',
        regex => 'fastp (\d+.\d+\.?\d+?)',
        desc  => 'Quality filtering of raw reads',
        paper => 'Chen S (2018)',
        doi   => 'https://doi.org/10.1093/bioinformatics/bty560',
        uri   => 'https://github.com/OpenGene/fastp',
    },
    'shovill' => {
        cmd   => 'shovill --version',
        regex => 'shovill (\d+.\d+\.?\d+?)',
        desc  => 'Assembly pf bacterial isolate genomes from Illumina paired-end reads',
        paper => '',
        uri   => 'https://github.com/tseemann/shovill',
    },   
    'mlst' => {
        cmd   => 'mlst --version',
        regex => 'mlst (\d+.\d+\.?\d+?)',
        desc  => 'Detection of MLST',
        paper => '',
        uri   => 'https://github.com/tseemann/mlst'
    },
    'unicycler' => {
        cmd   => 'unicycler --version',
        regex => 'Unicycler v(\d+.\d+\.?\d+?)',
        desc  => 'Hybrid assembly of short and long reads',
        paper => 'Wick RR (2017)',
        doi   => 'https://doi.org/10.1371/journal.pcbi.1005595',
        uri   => 'https://github.com/rrwick/Unicycler',
    },
    'kraken2' => {
        cmd   => 'kraken2 --version',
        regex => 'Kraken version (\d+.\d+\.?\d+?)',
        desc  => 'Taxonomic classification of sequencing reads',
        paper => 'Wood DE (2019)',
        doi   => 'https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0',
        uri   => 'https://github.com/DerrickWood/kraken2',
    },
    'abricate' => {
        cmd   => 'abricate --version',
        regex => 'abricate (\d+.\d+.?\.?\d+?)',
        desc  => 'Detection of virulence and resistome',
        paper => 'Seemann T, Abricate, Github https://github.com/tseemann/abricate',
        uri   => 'https://github.com/tseemann/abricate',
    },
    'multiqc' => {
        cmd   => 'multiqc --version',
        regex => 'multiqc, version (\d+.\d+)',
        desc  => 'Report generation',
        paper => 'Ewels P. et al., MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (2015)',
        uri   => 'https://multiqc.info/',
    },
    'snp-dists' => {
        cmd   => 'snp-dists -v',
        regex => '(\d+.\d+)',
        desc  => 'FASTA alignment to SNP distance matrix',
        uri   => 'https://github.com/tseemann/snp-dists',
    },
    'quast' => {
        cmd   => 'quast --version',
        regex => 'QUAST v(\d+.\d+\.?\d+?)',
        desc  => 'Quality assessment tool for genome assemblies',
        paper => 'Gurevich A (2013)',
        doi   => 'https://pubmed.ncbi.nlm.nih.gov/23422339/',
        uri   => 'http://bioinf.spbau.ru/quast',
    },
    
};

say ${HTML} qq(
<h1>GMH Nullarflow</h1>
<p class="lead">List of used software and links to their homepage or publication</p>
<dl class=dl-horizontal>);


for my $tool_name ( sort keys %{ $software } ) {
    my $tool = $software->{$tool_name};
    my $version = run($tool->{cmd}, $tool->{regex});
    say "$tool_name: $version";

    my $paper = '';
    my $desc  = $tool->{desc} // '...';
    if ($tool->{doi}) {
        my $link = $tool->{paper} // 'link';
        $paper = qq(&mdash; <a href="$tool->{doi}">$link</a>);
    }
    $tool_name = qq(<a href="$tool->{uri}">$tool_name</a>) if ($tool->{uri});

    say ${HTML}  qq(<dt>
    <strong>$tool_name</strong> version $version</dt>
    <dd>$desc  $paper
    </dd>);
 
}

sub run {
    my ($cmd, $regex, $canfail) = @_;
    my $out = `$cmd 2>&1`;
    if ( $? and  $canfail != 1 ) {
        die "Error executing '$cmd'\n";
    }

    if ($out =~/$regex/) {
        return $1;
    }
    return "unknown";
}