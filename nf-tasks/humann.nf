

nextflow.enable.dsl = 2
params.reads = "$baseDir/reads/*_R{1,2}.fastq.gz"
params.outdir = "nf-humann"

params.uniref = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/uniref/"
params.chocophlan = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/chocophlan/"
params.metaphlandb = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/mpa/"

reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

log.info """
         GMH Humann (version 0.1)
         ===================================
         input reads  : ${params.reads}
         outdir       : ${params.outdir}
         
         metaphlan    : ${params.metaphlandb}
         chocophlan   : ${params.chocophlan}
         uniref       : ${params.uniref}

         Running on   : ${params.max_memory}, ${params.max_cpus} cores
         """
         .stripIndent()
         
def UNIREF = file(params.uniref, checkIfExists: true)
def METAPHLANDB = file(params.metaphlandb, checkIfExists: true)
def CHOCOPHLAN = file(params.chocophlan, checkIfExists: true)

process VERSIONS {
    label "process_low"
    
    output:
    file("versions.txt")

    script:
    """
    fastp --version 2>&1 > versions.txt
    seqfu version 2>&1 >> versions.txt
    humann --version 2>&1 >> versions.txt
    """
}
process FASTP {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "filt $sample_id"
    label "process_low"

    input:
    tuple val(sample_id), path(reads) 
    path("versions.lock")
    
    output:
    tuple val(sample_id), path("${sample_id}.fq.gz")
  
    script:
    """
    mkdir -p filt
    fastp -i ${reads[0]} -I ${reads[1]} --stdout \
      --detect_adapter_for_pe --length_required 75 --thread ${task.cpus} \
      | seqfu cat --strip-name --strip-comments  --prefix "${sample_id}." \
      | gzip -c > ${sample_id}.fq.gz
    """     
    stub:
    """
    cat ${reads[0]} | head -n 4000 | gzip -c > ${sample_id}.fq.gz
    """
}
process INTERLEAVE {
    tag "ilv $sample_id"
    label "process_low"
    
    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}.fastq")
  
    script:
    """
    seqfu interleave -1 ${reads[0]} -2 ${reads[1]} | seqfu cat --strip-name --strip-comments  --prefix "${sample_id}." > ${sample_id}.fastq
    """  
    stub:
    """
    gzip -dc ${reads[0]} > ${sample_id}.fastq
    """
}  
process HUMANN {
    tag "$sample_id"
    publishDir "$params.outdir/humann/", mode:'copy'
    label "humann"
    
    input:
    tuple val(sample_id), path(reads) 
    path(chocophlan)
    path(uniref)
    path(metaphlandb)
    
    output:
    tuple val(sample_id), path("${sample_id}_genefamilies.tsv"), emit: genefamilies
    tuple val(sample_id), path("${sample_id}_pathabundance.tsv"), emit: pathabundance
    tuple val(sample_id), path("${sample_id}_pathcoverage.tsv"), emit: pathcoverage
    tuple val(sample_id), path("${sample_id}_metaphlan_bugs_list.tsv"), emit: metaphlan
    
  
    script:
    """
    INDEX=\$(cat mpa/mpa_latest)
    humann -i ${reads} -o $sample_id --threads ${task.cpus} --output-basename $sample_id \
        --nucleotide-database ${chocophlan} \
        --protein-database ${uniref} \
        --metaphlan-options="-x \$INDEX --bowtie2db ${metaphlandb} --unknown_estimation"
    mv ${sample_id}/*tsv .
    mv ${sample_id}/${sample_id}_humann_temp/${sample_id}_metaphlan_bugs_list.tsv ${sample_id}_metaphlan_bugs_list.tsv
    """  
    stub:
    """
    getoutput.py $sample_id
    """
}  
 
process multiqc {
    publishDir params.outdir, mode:'copy'
       
    input:
    path '*'  
    
    output:
    path 'multiqc_*'
     
    script:
    """
    multiqc . 
    """
} 

process JOIN {
    publishDir "$params.outdir", mode:'copy'
    input:
    path '*'

    output:
    path 'tables/*tsv'

    script:
    """
    mkdir -p tables
    humann_join_tables -i ./ --file_name _genefamilies.tsv -s -o tables/GeneFamilies.tsv
    humann_join_tables -i ./ --file_name _pathabundance.tsv -s -o tables/PathAbundance.tsv
    humann_join_tables -i ./ --file_name _pathcoverage.tsv -s -o tables/PathCoverage.tsv
    """
}

process JOINMPA {
    publishDir "$params.outdir", mode:'copy'
    input:
    path "*"

    output:
    path "tables/Taxa_merged.tsv"

    script:
    """
    mkdir -p tables
    merge_metaphlan_tables.py *_metaphlan_bugs_list.tsv > tables/Taxa_merged.tsv
    """
}

workflow {
   VERSIONS()
   FASTP(reads, VERSIONS.out)
   HUMANN(FASTP.out, CHOCOPHLAN, UNIREF, METAPHLANDB)
   JOINMPA((HUMANN.out.metaphlan).map{it -> it[1]}.collect())
   JOIN( HUMANN.out.genefamilies.mix( HUMANN.out.pathabundance, HUMANN.out.pathcoverage).map{it -> it[1]}.collect())
}