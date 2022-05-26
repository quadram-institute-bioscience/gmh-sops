

nextflow.enable.dsl = 2
params.reads = "$baseDir/reads/*_R{1,2}.fastq.gz"
params.outdir = "nf-humann"

params.uniref = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/uniref/"
params.chocophlan = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/chocophlan/"
params.metaphlandb = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/mpa/"

reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

log.info """
         GMH Humann (version 0.2)
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
    tuple val(sample_id), path("${sample_id}.fq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    
  
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} --stdout \
      --detect_adapter_for_pe --length_required 75 --thread ${task.cpus} \
      -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json \
      | seqfu cat --strip-name --strip-comments  --prefix "${sample_id}." \
      | gzip -c > ${sample_id}.fq.gz
    """     
    stub:
    """
    cat ${reads[0]} | head -n 4000 | gzip -c > ${sample_id}.fq.gz
    echo '{ "summary": { "after_filtering": { "total_reads":1090016, "total_bases":159862544,"read1_mean_length":146}}}' > ${sample_id}.fastp.json
    touch ${sample_id}.fastp.html
    """
}

process STATS {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "$sample_id"
    label "process_low"

    input:
    tuple val(sample_id), path(reads) 
    
    
    output:
    tuple val(sample_id), path("${sample_id}.tsv")
  
    script:
    """
    seqfu stats ${reads[0]} > ${sample_id}.tsv
    """     
    stub:
    """
    echo 'File,#Seq,Total bp,Avg,N50,N75,N90,auN,Min,Max' > tmp
    echo -n '${sample_id}' >> tmp
    echo ',16434,2465100,150.0,150,150,150,0.009,150,150' >> tmp
    sed 's/,/\\t/g' tmp > ${sample_id}.tsv
    """
}
process JOIN_STATS {
    input:
    path "*"

    output:
    path "stats.tsv"

    script:
    """
    cat *.tsv | head -n 1 > stats.tsv
    grep -v 'File' *.tsv | cut -f 2 -d ":" | sort >> stats.tsv
    """
}

process JSON_STATS {
    input:
    path "*"

    output:
    path "stats.tsv"

    script:
    """
    getCountsFromFastp.py --fields 10 --header --output stats.tsv *.fastp.json
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
    path 'tables/*.tsv'

    script:
    """
    mkdir -p tables
    humann_join_tables -i ./ --file_name _genefamilies.tsv -s -o tables/GeneFamilies.tsv
    humann_join_tables -i ./ --file_name _pathabundance.tsv -s -o tables/PathAbundance.tsv
    humann_join_tables -i ./ --file_name _pathcoverage.tsv -s -o tables/PathCoverage.tsv
    """
}



process TOP_TAXA {
    label 'process_filtering'
    publishDir "$params.outdir/top_taxa/", mode:'copy'

    input:
    path "merged-taxonomy.tsv"

    output:
    path '*.csv'

    script:
    """
    #Estimate the top N genefamilies and pathways in the study
    extract_taxaabundance.py -i merged-taxonomy.tsv -o Species -r S -t 20
    extract_taxaabundance.py -i merged-taxonomy.tsv -o Genus -r G -t 20
    """
}
process JOIN_TAXONOMY {
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

process TOP_PATHWAYS {
    publishDir "$params.outdir/", mode:'copy'
    input:
    path tables
    path 'stats.tsv'

    output:
    path 'summary'
    
    script:
    """
    #Estimate the top N genefamilies and pathways in the study
    mkdir -p summary
    extract_genefamilies.py -i GeneFamilies.tsv -st stats.tsv -qc summary/qc -o summary/top20 -t 20
    extract_pathabundance.py -i PathAbundance.tsv -o summary/top20 -t 20
    """
}

process BUBBLE_PLOTS {
    publishDir "$params.outdir/", mode:'copy'

    input:
    path summary

    output:
    path 'plots-pathways'

    script:
    """
    #Plot the top N genefamilies and pathways in the study
    bubble_plot.R summary/qc_report.csv summary/top20_GFabundance-rel.csv top20gf False ./
    bubble_plot.R summary/qc_report.csv summary/top20_PWabundance-rel.csv top20pa False ./
    mv plots plots-pathways
    """
}
process BUBBLE_TAXA {
    publishDir "$params.outdir/", mode:'copy'
    input:
    path top_taxa

    output:
    path 'plots-taxa'

    script:
    """
    #Plot the top N genefamilies and pathways in the study
    mkdir -p plots
    heatmap.R Species-rel.csv ./plots-taxa/
    heatmap.R Genus-rel.csv ./plots-taxa
    """
}
workflow {
   // Collect versions (and ensures that the programs are available)
   VERSIONS()
   // Clean reads, interleave (versions is used as lock to run VERSION first)
   FASTP(reads, VERSIONS.out)

   // Count reads, merge counts
   //STATS(FASTP.out.reads)
   //JOIN_STATS( STATS.out.map{it -> it[1]}.collect() )
   JSON_STATS(FASTP.out.json.map{it -> it[1]}.collect())

   // RUN HUMANN and use its tools to merge tables
   HUMANN(FASTP.out.reads, CHOCOPHLAN, UNIREF, METAPHLANDB)
   JOIN_TAXONOMY((HUMANN.out.metaphlan).map{it -> it[1]}.collect())
   JOIN( HUMANN.out.genefamilies.mix( HUMANN.out.pathabundance, HUMANN.out.pathcoverage).map{it -> it[1]}.collect())

   // Sumeet Tiwari scripts to summarise results
   TOP_PATHWAYS(JOIN.out, JSON_STATS.out)
   BUBBLE_PLOTS(TOP_PATHWAYS.out)
   TOP_TAXA(JOIN_TAXONOMY.out)
   BUBBLE_TAXA(TOP_TAXA.out)
}
