#!/usr/bin/env nextflow

version = '0.1.3__2021b'

def helpMessage() {

    log.info"""
    Usage:
    nextflow run magenta.nf --reads '*_R{1,2}.fastq.gz' 

    Main arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --outdir                      Output directory (default ./metagenome)
      --tempdir                     Absolute PATH to the temporary directory (default ./tmp)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, singularity, test and more.

    Databases:
      --kraken2db
      --kraken2secondary
      --metaphlandb
      --chocophlan
      --uniref
    
    Skipping pre-enabled steps:
      --skip_taxa    
        --skip_kraken
        --skip_kraken_secondary
        --skip_metaphlan
      --skip_denovo
        --skip_binning
        --skip_eggs
    
    Include optional steps:
      --virome

    """.stripIndent()
}
params.skip_kraken = false
params.skip_kraken_secondary = false
params.skip_eggs = false
params.skip_binning = false
params.skip_taxa = false
params.skip_denovo = false

params.readPaths = false
params.singleEnd = false
params.reads = ""
params.cpus = 8
params.skip_eggs = false
params.outdir = "assemblies"
params.tempdir = "/tmp/"
params.adapter_forward  = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse  = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality     = 20
params.trimming_quality = 14
params.min_length       = 80
params.subsample        = 0.6

params.metaphlandb      = '/qib/platforms/Informatics/transfer/outgoing/gmh/databases/metaphlan_databases/'
params.chocophlan       = '/qib/platforms/Informatics/transfer/outgoing/gmh/databases/chocophlan/'
params.uniref           = '/qib/platforms/Informatics/transfer/outgoing/gmh/databases/uniref/'

params.kraken2db        = '/qib/platforms/Informatics/transfer/outgoing/databases/kraken2/mini'
params.kraken2secondary = ''

params.virome = false

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


if(!params.skip_kraken){
    Channel
        .fromPath( "${params.kraken2db}", checkIfExists: true )
        .set { file_kraken2db }     
} else {
    file_kraken2db = Channel.from()
}

if ( !params.skip_kraken &&  !params.skip_kraken_secondary ) {
    Channel
        .fromPath( "${params.kraken2secondary}", checkIfExists: true )
        .set { file_kraken2db_secondary }     
} else {
    file_kraken2db_secondary = Channel.from()
}


Channel
        .fromPath( "${params.metaphlandb}", checkIfExists: true )
        .set { file_metaphlandb }    

Channel
        .fromPath( "${params.chocophlan}", checkIfExists: true )
        .set { file_chocophlan }    

Channel
        .fromPath( "${params.uniref}", checkIfExists: true )
        .set { file_uniref }  

log.info """
===        GMH - VirOne v.${version}          ===
============================================
 reads        : ${params.reads}
 tempdir      : ${params.tempdir}
 outdir       : ${params.outdir}
 kraken2 A    : ${params.kraken2db}
 kraken2 B    : ${params.kraken2secondary}
 Metaphlan    : ${params.metaphlandb}
                ${params.chocophlan}
                ${params.uniref}
 """


process versions {
   publishDir params.outdir, mode: "copy"
   label 'onecore'
   executor 'local'

   output:
   file('versions.html') into versions_ch

   script:
   """
   checker.py --html versions.html --filepath "${params.uniref}" "${params.tempdir}" "${params.metaphlandb}" "${params.chocophlan}"

   if [ ! -e "${params.metaphlandb}/mpa_latest" ];
   then
      if [ -e "${params.metaphlandb}/mpa_previous" ]; then
        mv "${params.metaphlandb}/mpa_previous" "${params.metaphlandb}/mpa_latest"
      else
        echo " ERROR: Unable to find metaphlan database"
      fi
   fi
   """
}


process filter {
    publishDir "${params.outdir}/reads/", mode: "copy"
    tag "$name"
    label 'lowmem'

    input:
    file('versions.txt') from versions_ch
    set val(name), file(reads)   from read_pairs_ch
    val trim_qual from params.trimming_quality
    val adapter from params.adapter_forward
    val adapter_reverse from params.adapter_reverse
    val qual from params.mean_quality
    val min_len from params.min_length

    output:
    set val(name), file("${name}.filtered*.fastq.gz") into cleaned_for_prokka, cleaned_for_kraken, cleaned_for_mapping, cleaned_for_metaphlan, cleaned_for_binning
    set val(name), file("*.fastp.*") into json_ch

    //TODO - Add kraken2 remove host?

    script:
    //def pe_input   = params.singleEnd ? '' :  "-I \"${reads[1]}\""
    //def pe_output1 = params.singleEnd ? "-o \"${name}.filtered.fastq.gz\"" :  "-o \"${name}.filtered_R1.fastq.gz\""
    //def pe_output2 = params.singleEnd ? '' :  "-O \"${name}.filtered_R2.fastq.gz\""
    if (params.singleEnd)
      """
      # Single End
      fastp -w "${task.cpus}" -q "${qual}"\
          --cut_by_quality3  --cut_by_quality5  \
          --cut_mean_quality "${trim_qual}" \
          --json ${name}.fastp.json \
          --adapter_sequence=${adapter}  \
          -i "${reads[0]}"   \
          -o "${name}.filtered.fastq.gz"

      """
    else
      """
      # Paired-End [was:  --adapter_sequence=${adapter} --adapter_sequence_r2=${adapter_reverse} ]

      fastp -w "${task.cpus}" --qualified_quality_phred "${qual}" \
          --cut_by_quality5 --cut_by_quality3  \
          --cut_mean_quality "${trim_qual}" \
          --length_required "${min_len}" \
          --json ${name}.fastp.json \
          --detect_adapter_for_pe \
          -i "${reads[0]}"  -I "${reads[1]}" \
          -o "${name}.filtered_R1.fastq.gz" -O "${name}.filtered_R2.fastq.gz" 
      """

}


process kraken2 {
    publishDir "${params.outdir}/taxonomy/", mode: "copy"
    tag "$name"

    when:
    !params.skip_kraken && !params.skip_taxa

    input:
    set val(name), file(reads)   from cleaned_for_kraken
    file(db1) from file_kraken2db
    file(db2) from file_kraken2db_secondary

    output:
    set val(name), file("*.report.txt") into kraken2_ch

    script:
    if ( params.skip_kraken_secondary  )
      """
      kraken2 --threads ${task.cpus} --db ${db1}/  --report ${name}.${db1.baseName}.report.txt --paired --confidence 0.5 "${reads[0]}" "${reads[1]}" > /dev/null
      """
    else
      """
      kraken2 --threads ${task.cpus} --db ${db1}/  --report ${name}.${db1.baseName}.report.txt --paired --confidence 0.5 "${reads[0]}" "${reads[1]}" > /dev/null
      kraken2 --threads ${task.cpus} --db ${db2}/  --report ${name}.${db2.baseName}.report.txt --paired --confidence 0.5 "${reads[0]}" "${reads[1]}" > /dev/null
      """

}

process metaphlan {
  
    publishDir "${params.outdir}/taxonomy/", mode: "copy"
    tag "$name"

    when:
    !params.skip_taxa
    
    input:
    set val(name), file(reads)   from cleaned_for_metaphlan
    file(metaphlandb) from file_metaphlandb
    file(uniref) from file_uniref
    file(chocophlan) from file_chocophlan

    output:
    set val(name), file("*.tsv") into metaphlan_ch

    script:
    """
    logjob | tee job.md
    cat ${reads[0]} ${reads[1]} | gzip -d > ${name}.fq

    # Run humann / metaphlan
    humann \
        -i ${name}.fq \
        -o out \
        --nucleotide-database ${chocophlan} \
        --protein-database ${uniref} \
        --metaphlan-options='--bowtie2db ${metaphlandb}  --unknown_estimation' \
        --threads ${task.cpus}
   

    humann_renorm_table --input out/${name}_genefamilies.tsv \
      --output ${name}_genefamilies-cpm.tsv --units cpm --update-snames
    humann_renorm_table --input out/${name}_genefamilies.tsv \
      --output ${name}_genefamilies-relab.tsv --units relab --update-snames
    humann_renorm_table --input out/${name}_pathabundance.tsv \
      --output ${name}_pathabundance-cpm.tsv --units cpm --update-snames
    humann_renorm_table --input out/${name}_pathabundance.tsv \
      --output ${name}_pathabundance-relab.tsv --units relab --update-snames

    mv out/*/*_metaphlan_bugs_list.tsv out/${name}_metaphlan.tsv
    """

}

process megahit {
    publishDir "${params.outdir}/assembly/", mode: "copy"

    tag "$name"

    when:
    !(params.skip_denovo)

    input:
    set val(name), file(reads)   from cleaned_for_prokka

    output:
    set val(name), file("${name}.fasta") into contigs_ch, contigs_virsorter1, contigs_map_ch, contigs_bin_ch, contigs_das_ch

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} -o ctg -t ${task.cpus} 
    { fu-rename --prefix ${name} --separator _ ctg/final.contigs.fa --nocomments | seqkit seq --only-id > ${name}.fasta; } || mv ctg/final.contigs.fa > ${name}.fasta
    """

}

process virfinder1 {
  publishDir "${params.outdir}/virome/", mode: "copy"

  tag "$name"

  input:
  set val(name), file('contigs.fasta')          from contigs_virsorter1

  when:
  !(params.skip_denovo) && params.virome
    
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

process prodigal {
    publishDir "${params.outdir}/annotation/", mode: "copy"

    tag "$name"

    when:
    !(params.skip_denovo)
    
    input:
    set val(name), file('contigs.fasta')   from contigs_ch

    output:
    set val(name), file("*.{faa,ffn}") into prokka_prot
    set val(name), file("*.gff") into prokka_gff


    script:
    """
    prodigal -a ${name}.faa -f gff -d ${name}.ffn -p meta -o ${name}.gff -i contigs.fasta
    """

}

process eggnog {
  publishDir "${params.outdir}/annotation/", mode: "copy"
  tag "$name"

  when:
  !(params.skip_denovo) && !(params.skip_eggs)

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
    publishDir "${params.outdir}/align/", mode: "copy"
    tag "$name"
    label 'hicpu'

    when:
    !(params.skip_denovo)
    
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
        rpkm=${name}.mapping.rpkm covhist=${name}.mapping.hist trimreaddescriptions=true

    samtools view -bS ${name}.unsorted.sam | samtools fixmate - fixed.bam
    samtools sort -@ ${task.cpus} -o ${name}.mapping.bam fixed.bam
    samtools index ${name}.mapping.bam

    rm ${name}.unsorted.sam fixed.bam
    """

}

process coverage {
  publishDir "${params.outdir}/align/", mode: "copy"
  tag "$name"
  label 'lowmem'

  when:
  !(params.skip_denovo)
    
  input:
  set val(name), file(bam)               from alncov_ch
  set val(name), file("ann.gff")         from prokka_gff

  output:
  set val(name), file("${name}.bed")          into cov_ch
  set val(name), file("${name}.counts")       into counts_ch
  script:
  """
  #mkdir -p ${params.tempdir}
  #featureCounts -T ${task.cpus} -Q 1 --tmpDir ${params.tempdir} -a ann.gff -o ${name}.counts -t CDS -g ID ${bam[0]}
  covtobed --discard-invalid-alignments --min-cov=1 --physical-coverage ${bam[0]} | tee ${name}.bed | covtotarget ann.gff > ${name}.counts

  """
}

process binning {
    publishDir "${params.outdir}/binning/", mode: "copy"

    tag "$name"
    
    when:
    !(params.skip_denovo)
    
    input:
    set val(name), file(bam)               from aln_ch
    set val(name), file(reads)             from cleaned_for_binning
    set val(name), file('contigs.fasta')   from contigs_bin_ch

    output:
    set val(name), file("*metabat*") into metabat_ch
    set val(name), file("*maxbin*") into maxbin_ch

    script:
    """
    echo BINNING
    echo ID  : ${name} - ${reads}
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
  publishDir "${params.outdir}/binning/", mode: "copy"

  when:
  !(params.skip_denovo)
    
  tag "$name"

  input:
  set val(name), file('contigs.fna')   from contigs_das_ch
  set val(name), file(metabat) from metabat_ch
  set val(name), file(maxbin)  from  maxbin_ch

  output:
  set val(name), file("*summary.txt"), file("*.pdf") optional true into dasqc_ch
  set val(name), file("refine_DASTool_bins/*.fa")    optional true into das_ch

  script:
  """
  BinningStep.pl ${name}
  """
}

workflow.onComplete {
    log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    log.info ( workflow.success ?
        "\nDone! The results are saved in --> $params.outdir/\n" :
        "Oops .. something went wrong" )
}
