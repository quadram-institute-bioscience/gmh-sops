

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
    file("versions.txt")
    
    output:
    tuple val(sample_id), path("filt/${sample_id}_R{1,2}.fastq.gz")
  
    script:
    """
    mkdir -p filt
    fastp -i ${reads[0]} -I ${reads[1]} -o filt/${sample_id}_R1.fastq.gz -O filt/${sample_id}_R2.fastq.gz \
      --detect_adapter_for_pe --length_required 75 --thread ${task.cpus}
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
}  
process HUMANN {
    tag "$sample_id"
    publishDir params.outdir, mode:'copy'
    label "humann"
    
    input:
    tuple val(sample_id), path(reads) 
    path(chocophlan)
    path(uniref)
    path(metaphlandb)
    
    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_genefamilies.tsv"), emit: genefamilies
    tuple val(sample_id), path("${sample_id}/${sample_id}_pathabundance.tsv"), emit: pathabundance
    tuple val(sample_id), path("${sample_id}/${sample_id}_pathcoverage.tsv"), emit: pathcoverage
    tuple val(sample_id), path("${sample_id}/${sample_id}_humann_temp/${sample_id}_metaphlan_bugs_list.tsv"), emit: metaphlan
    
  
    script:
    """
    humann -i ${reads} -o $sample_id --threads ${task.cpus} --output-basename $sample_id \
        --nucleotide-database ${chocophlan} \
        --protein-database ${uniref} \
        --metaphlan-options="--bowtie2db ${metaphlandb} --add_viruses --unknown_estimation"

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
    multiqc --cl_config "prokka_fn_snames: True" . 
    """
} 

workflow {
   VERSIONS()
   FASTP(reads, VERSIONS.out)
   INTERLEAVE(FASTP.out)
   HUMANN(INTERLEAVE.out, CHOCOPHLAN, UNIREF, METAPHLANDB)
    
}