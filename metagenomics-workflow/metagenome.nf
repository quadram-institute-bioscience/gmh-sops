#!/usr/bin/env nextflow

version = '0.1.1'

def helpMessage() {

    log.info"""
    Usage:
    nextflow run metagmh --reads '*_R{1,2}.fastq.gz' -profile docker

    Main arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --outdir                      Output directory (default ./metagenome)
      --tempdir                     Absolute PATH to the temporary directory (default ./tmp)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, singularity, test and more.

    Databases:
      --kraken2db
      --kraken2secondary
    Skipping steps:
      --skip-kraken
      --skip-metaphlan
      --skip-binning
      --skip-eggs

    """.stripIndent()
}

params.readPaths = false
params.reads = ""
params.cpus = 8
params.skipeggs = false
params.outdir = "assemblies"
params.tempdir = "/tmp/"
params.adapter_forward  = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse  = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality     = 18
params.trimming_quality = 14

params.subsample        = 0.6
params.kraken2db        = '/qi/db/kraken2/kraken2_db_20190111'

if (params.readPaths) {
    if (params.singleEnd) {

        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_pairs_ch; read_pairs2_ch;  }
    } else {

        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_pairs_ch; read_pairs2_ch;  }
    }
} else {

    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { read_pairs_ch; read_pairs2_ch;  }
}




log.info """
===        GMH - VirOne v.${version}          ===
============================================
 reads        : ${params.reads}
 tempdir      : ${params.tempdir}
 outdir       : ${params.outdir}
 kraken2 A    : ${params.kraken2db}
 kraken2 B    : ${params.kraken2secondary}
 """

 process versions {
   publishDir params.outdir, mode: "copy"
   label 'onecore'

   output:
   file('versions.txt') into versions_ch

   script:
   """
   date > versions.txt
   echo '# FASTP'       >> versions.txt
   fastp -v 2>&1        >> versions.txt
   echo '# SAMTOOLS'    >> versions.txt
   samtools --version   >> versions.txt
   echo '# MEGAHIT'     >> versions.txt
   megahit --version    >> versions.txt
   echo '# MAXBIN'      >> versions.txt
   run_MaxBin.pl  -version  >> versions.txt
   echo '# METABAT'     >> versions.txt
   metabat2 |& grep version >> versions.txt || true
   echo '# SEQFU'       >> versions.txt
   fu-cat --version     >> versions.txt
   echo '# KRAKEN2'     >> versions.txt
   kraken2 --version |& grep Kraken >> versions.txt
   zip
   """
 }

process filter {
    publishDir params.outdir, mode: "copy"
    tag "$name"
    label 'lowmem'

    input:
    file('versions.txt') from versions_ch
    set val(name), file(reads)   from read_pairs_ch
    val trim_qual from params.trimming_quality
    val adapter from params.adapter_forward
    val adapter_reverse from params.adapter_reverse
    val qual from params.mean_quality


    output:
    set val(name), file("${name}.filtered*.fastq.gz") into cleaned_for_prokka, cleaned_for_kraken, cleaned_for_mapping, cleaned_for_metaphlan, cleaned_for_binning
    set val(name), file("*.fastp.*") into json_ch

    script:
    def pe_input   = params.singleEnd ? '' :  "-I \"${reads[1]}\""
    def pe_output1 = params.singleEnd ? "-o \"${name}.filtered.fastq.gz\"" :  "-o \"${name}.filtered_R1.fastq.gz\""
    def pe_output2 = params.singleEnd ? '' :  "-O \"${name}.filtered_R2.fastq.gz\""
    """
    mkdir -p ${params.tempdir}
    logjob | tee job.md

    fastp -w "${task.cpus}" -q "${qual}" --cut_by_quality5 \
        --cut_by_quality3 --cut_mean_quality "${trim_qual}" \
        --json ${name}.fastp.json \
        --adapter_sequence=${adapter} --adapter_sequence_r2=${adapter_reverse} \
        -i "${reads[0]}" $pe_input $pe_output1 $pe_output2

    """
}


process kraken2 {
    publishDir params.outdir, mode: "copy"
    tag "$name"


    input:
    set val(name), file(reads)   from cleaned_for_kraken

    output:
    set val(name), file("*_report.txt") into kraken2_ch

    script:
    """
    logjob | tee job.md
    kraken2 --threads ${task.cpus} --db ${params.kraken2db}         --report ${name}.kraken2_report.txt --paired --confidence 0.5 "${reads[0]}" "${reads[1]}" > /dev/null
    kraken2 --threads ${task.cpus} --db ${params.kraken2secondary}  --report ${name}.gtdbk_report.txt --paired --confidence 0.5 "${reads[0]}" "${reads[1]}" > /dev/null
    """

}

process metaphlan {
    publishDir params.outdir, mode: "copy"
    tag "$name"

    input:
    set val(name), file(reads)   from cleaned_for_metaphlan

    output:
    set val(name), file("*.metaphlan3*") into metaphlan_ch

    script:
    """
    logjob | tee job.md
    metaphlan --input_type fastq ${reads[0]},${reads[1]} \
                --bowtie2out bowtie_aln.bz2 \
                --bowtie2db ${params.metaphlandb} \
                --nproc ${task.cpus} \
                --add_viruses \
                --unknown_estimation \
                --biom ${name}.metaphlan3.biom  \
                -o ${name}.metaphlan3.txt

    # Fail safe
    touch  ${name}.metaphlan3.biom  ${name}.metaphlan3.txt
    """

}

process megahit {
    publishDir params.outdir, mode: "copy"

    tag "$name"

    input:
    set val(name), file(reads)   from cleaned_for_prokka

    output:
    set val(name), file("${name}.fasta") into contigs_ch, contigs_virsorter1, contigs_map_ch, contigs_bin_ch, contigs_das_ch

    script:
    """
    if [[ ! -d "${params.tempdir}" ]];
    then
      mkdir -p ${params.tempdir};
    fi
    logjob | tee job.md
    megahit -1 ${reads[0]} -2 ${reads[1]} -o ctg -t ${task.cpus} --tmp-dir ${params.tempdir}
    { fu-rename --prefix ${name} --separator _ ctg/final.contigs.fa --nocomments | seqkit seq -i > ${name}.fasta; } || mv ctg/final.contigs.fa > ${name}.fasta
    """

}

process virfinder1 {
  publishDir params.outdir, mode: "copy"

  tag "$name"

  input:
  set val(name), file('contigs.fasta')          from contigs_virsorter1

  output:
  set val(name), file("${name}.virsorter1.csv") into virsorter1_csv
  set val(name), file("${name}.virsorter1.zip") into virsorter1_zip

  script:
  """
  # db=1 (Refseqdb), db=2 (Viromedb)
  wrapper_phage_contigs_sorter_iPlant.pl --fna contigs.fasta --wdir \$PWD/virfinder1  --db 1 --ncpu ${task.cpus}
  mv virfinder1/VIRSorter_global-phage-signal.csv ${name}.virsorter1.csv
  mv virfinder1/Predicted_viral_sequences ${name}.predicted_viral_sequences
  zip -r ${name}.virsorter1.zip ${name}.predicted_viral_sequences
  """
}

process metaprokka {
    publishDir params.outdir, mode: "copy"

    tag "$name"

    input:
    set val(name), file('contigs.fasta')   from contigs_ch

    output:
    set val(name), file("*.{faa,ffn}") into prokka_prot
    set val(name), file("*.log") into prokka_log
    set val(name), file("*.gff") into prokka_gff


    script:
    """
    logjob | tee job.md
    metaprokka --metagenome --cpu ${task.cpus} --locustag ${name} --outdir prokka --prefix annotation contigs.fasta
    mv prokka/annotation.gff ${name}.gff
    mv prokka/annotation.faa ${name}.faa
    mv prokka/annotation.ffn ${name}.ffn
    mv prokka/annotation.log ${name}.prokka.log
    """

}

process eggnog {
  publishDir params.outdir, mode: "copy"
  tag "$name"

  when:
  !params.skipeggs

  input:
  set val(name), file(prokka)            from prokka_prot

  output:
  set val(name), file("*.annotations.tsv")   into eggnog_ch

  script:
  """
  logjob | tee job.md
  #{out}.emapper.seed_orthologs {out}.emapper.annotations
  emapper.py -i "${prokka[0]}" --output ${name} -m diamond --cpu ${task.cpus}
  rm "${name}.emapper.seed_orthologs"
  mv "${name}.emapper.annotations" "${name}.annotations.tsv"
  """
}
process align {
    publishDir params.outdir, mode: "copy"
    tag "$name"
    label 'hicpu'

    input:
    set val(name), file(reads)             from cleaned_for_mapping
    set val(name), file('contigs.fasta')   from contigs_map_ch

    output:
    set val(name), file("*.mapping.bam*")  into aln_ch, alncov_ch
    set val(name), file("*hist"), file("*.rpkm") into alnstats_ch


    script:
    """
    bbmap.sh  in=${reads[0]} in2=${reads[1]} ref=contigs.fasta nodisk \
        threads=${task.cpus} minid=0.90 mappedonly=t out=${name}.unsorted.sam  \
        rpkm=${name}.mapping.rpkm covhist=${name}.mapping.hist

    samtools view -bS ${name}.unsorted.sam | samtools sort -@ ${task.cpus} -o ${name}.mapping.bam -
    samtools index ${name}.mapping.bam

    rm ${name}.unsorted.sam
    """

}

process coverage {
  publishDir params.outdir, mode: "copy"
  tag "$name"
  label 'lowmem'

  input:
  set val(name), file(bam)               from alncov_ch
  set val(name), file("ann.gff")         from prokka_gff

  output:
  set val(name), file("${name}.counts")          into cov_ch

  script:
  """
  mkdir -p ${params.tempdir}
  featureCounts -T ${task.cpus} -Q 1 --tmpDir ${params.tempdir} -a ann.gff -o ${name}.counts -t CDS -g ID ${bam[0]}

  covtobed --discard-invalid-alignments --min-cov=1  ${bam[0]} > coverage.bed

  """
}

process binning {
    publishDir params.outdir, mode: "copy"

    tag "$name"

    input:
    set val(name), file(bam)               from aln_ch
    set val(name), file(reads)             from cleaned_for_binning
    set val(name), file('contigs.fasta')   from contigs_bin_ch

    output:
    set val(name), file("*metabat*") into metabat_ch
    set val(name), file("*maxbin*") into maxbin_ch

    script:
    """
    echo BAM0: ${bam[0]}
    echo BAM1: ${bam[1]}

    jgi_summarize_bam_contig_depths \
        --outputDepth ${name}.depth.txt \
        ${name}.mapping.bam

    metabat2 \
        -i contigs.fasta \
        -a ${name}.depth.txt \
        -o ${name}.metabat \
        -t ${task.cpus}  -m 1500  -v --unbinned

    run_MaxBin.pl  \
        -contig contigs.fasta \
        -out ./${name}.maxbin \
        -thread ${task.cpus} \
        -reads ${reads[0]} \
        -reads2 ${reads[1]}
    """
}

process dastool {
  publishDir params.outdir, mode: "copy"

  tag "$name"

  input:
  set val(name), file('contigs.fna')   from contigs_das_ch
  set val(name), file(metabat) from metabat_ch
  set val(name), file(maxbin)  from  maxbin_ch

  output:
  set val(name), file("*.summary.txt"), file("*.pdf") into dasqc_ch
  set val(name), file("*.das.bin*.fa") into das_ch

  script:
  """
  # METABAT
  rm *unbinned.fa *tooShort.fa *lowDepth.fa

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
  """
}

workflow.onComplete {
    log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    log.info ( workflow.success ?
        "\nDone! The results are saved in --> $params.outdir/\n" :
        "Oops .. something went wrong" )
}
