#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mag
========================================================================================

nf-core/mag Analysis Pipeline. Started 2018-05-22.
#### Homepage / Documentation
https://github.com/nf-core/mag
#### Authors
Hadrien Gourlé HadrienG <hadrien.gourle@slu.se> - hadriengourle.com>
Daniel Straub <d4straub@gmail.com>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/nullarbor --reads '*_R{1,2}.fastq.gz' -profile docker
    nextflow run nf-core/nullarbor --manifest 'manifest.tsv' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/mltiqc.yaml"
params.email = false
params.plaintext_email = false
params.manifest = false
params.busco_reference = "https://busco-archive.ezlab.org/v3/datasets/bacteria_odb9.tar.gz"

ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}
// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * short read preprocessing options
 */
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false
// params.phix_reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Enterobacteria_phage_phiX174_sensu_lato/all_assembly_versions/GCA_002596845.1_ASM259684v1/GCA_002596845.1_ASM259684v1_genomic.fna.gz"
params.phix_reference = "$baseDir/assets/data/GCA_002596845.1_ASM259684v1_genomic.fna.gz"

 
/*
 * assembly options
 */
params.skip_spades = false
params.skip_spadeshybrid = false
params.skip_megahit = false
params.skip_quast = false

/*
 * taxonomy options
 */
params.kraken2_db = false


/*
 * Create a channel for input read files
 */
if(!params.skip_busco){
    Channel
        .fromPath( "${params.busco_reference}", checkIfExists: true )
        .set { file_busco_db }
} else {
    file_busco_db = Channel.from()
}

if(params.centrifuge_db){
    Channel
        .fromPath( "${params.centrifuge_db}", checkIfExists: true )
        .set { file_centrifuge_db }
} else {
    file_centrifuge_db = Channel.from()
}

if(params.kraken2_db){
    Channel
        .fromPath( "${params.kraken2_db}", checkIfExists: true )
        .set { file_kraken2_db }
} else {
    file_kraken2_db = Channel.from()
}

if(params.cat_db){
    Channel
        .fromPath( "${params.cat_db}", checkIfExists: true )
        .set { file_cat_db }
} else {
    file_cat_db = Channel.from()
}

if(!params.keep_phix) {
    Channel
        .fromPath( "${params.phix_reference}", checkIfExists: true )
        .set { file_phix_db }
}

def returnFile(it) {
// Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "Missing file in TSV file: ${inputFile}, see --help for more information"
    return inputFile
}

if(params.manifest){
    manifestFile = file(params.manifest)
    // extracts read files from TSV and distribute into channels
    Channel
        .from(manifestFile)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .splitCsv(sep:'\t')
        .map { row ->
            def id = row[0]
            def lr = returnFile( row[1] )
            def sr1 = returnFile( row[2] )
            def sr2 = returnFile( row[3] )
            [ id, lr, sr1, sr2 ]
            }
        .into { files_sr; files_all_raw }
    // prepare input for fastqc
    files_sr
        .map { id, lr, sr1, sr2 -> [ id, [ sr1, sr2 ] ] }
        .into { read_files_fastqc; read_files_fastp }

} else if(params.readPaths){
     if(params.singleEnd){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_fastp }
        files_all_raw = Channel.from()
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_fastp }
        files_all_raw = Channel.from()
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_fastp }
    files_all_raw = Channel.from()
 }



// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']                = params.reads
summary['Fasta Ref']            = params.fasta
summary['Data Type']            = params.singleEnd ? 'Single-End' : 'Paired-End'
if(params.centrifuge_db) summary['Centrifuge Db']   = params.centrifuge_db
if(params.kraken2_db) summary['Kraken2 Db']         = params.kraken2_db
if(!params.skip_busco) summary['Busco Reference']   = params.busco_reference
summary['Max Resources']        = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']           = params.outdir
summary['Launch dir']           = workflow.launchDir
summary['Working dir']          = workflow.workDir
summary['Script dir']           = workflow.projectDir
summary['User']                 = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']        = params.awsregion
   summary['AWS Queue']         = params.awsqueue
}
summary['Config Profile']       = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']     = params.email
  summary['MultiQC maxsize']    = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-mag-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/mag Workflow Summary'
    section_href: 'https://github.com/nf-core/mag'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.manifest.nextflowVersion > v_nextflow.txt
    fastp -v 2> v_fastp.txt
    kraken2 -v > v_kraken2.txt
    quast -v > v_quast.txt

    #scrape_software_versions.py > software_versions_mqc.yaml
    touch software_versions_mqc.yaml
    """
}


/*
 * Trim adapter sequences on long read nanopore files
 */
if (!params.skip_adapter_trimming) {
    process porechop {
        tag "$id"

        input:
        set id, file(lr), sr1, sr2 from files_all_raw

        output:
        set id, file("${id}_porechop.fastq"), sr1, sr2 into files_porechop
        set id, file(lr), val("raw") into files_nanoplot_raw

        script:
        """
        porechop -i ${lr} -t "${task.cpus}" -o ${id}_porechop.fastq
        """
    }
} else {
    files_all_raw
        .into{ files_porechop; pre_files_nanoplot_raw }
    pre_files_nanoplot_raw
        .map { id, lr, sr1, sr2 -> [ id, lr, "raw" ] }
        .set { files_nanoplot_raw }
}

/*
 * Remove reads mapping to the lambda genome.
 * TODO: add lambda phage to igenomes.config?
 */
if (!params.keep_lambda) {
    Channel
        .fromPath( "${params.lambda_reference}", checkIfExists: true )
        .set { file_nanolyse_db }
    process nanolyse {
        tag "$id"

        publishDir "${params.outdir}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "QC_longreads/NanoLyse/$filename" : null}

        input:
        set id, file(lr), file(sr1), file(sr2), file(nanolyse_db) from files_porechop.combine(file_nanolyse_db)

        output:
        set id, file("${id}_nanolyse.fastq.gz"), file(sr1), file(sr2) into files_nanolyse
        file("${id}_nanolyse_log.txt")

        script:
        """
        cat ${lr} | NanoLyse --reference $nanolyse_db | gzip > ${id}_nanolyse.fastq.gz

        echo "NanoLyse reference: $params.lambda_reference" >${id}_nanolyse_log.txt
        cat ${lr} | echo "total reads before NanoLyse: \$((`wc -l`/4))" >>${id}_nanolyse_log.txt
        zcat ${id}_nanolyse.fastq.gz | echo "total reads after NanoLyse: \$((`wc -l`/4))" >>${id}_nanolyse_log.txt
        """
    }
} else {
    files_porechop
        .set{ files_nanolyse }
}


/*
 * Quality filter long reads focus on length instead of quality to improve assembly size
 */
process filtlong {
    tag "$id"

    input:
    set id, file(lr), file(sr1), file(sr2) from files_nanolyse

    output:
    set id, file("${id}_lr_filtlong.fastq.gz") into files_lr_filtered
    set id, file("${id}_lr_filtlong.fastq.gz"), val('filtered') into files_nanoplot_filtered

    script:
    """
    filtlong \
        -1 ${sr1} \
        -2 ${sr2} \
        --min_length ${params.longreads_min_length} \
        --keep_percent ${params.longreads_keep_percent} \
        --trim \
        --length_weight ${params.longreads_length_weight} \
        ${lr} | gzip > ${id}_lr_filtlong.fastq.gz
    """
}


/*
 * Quality check for nanopore reads and Quality/Length Plots
 */
process nanoplot {
    tag "$id"
    publishDir "${params.outdir}/QC_longreads/NanoPlot_${id}", mode: 'copy'

    input:
    set id, file(lr), type from files_nanoplot_raw.mix(files_nanoplot_filtered)

    output:
    file '*.png'
    file '*.html'
    file '*.txt'

    script:
    """
    NanoPlot -t "${task.cpus}" -p ${type}_  --title ${id}_${type} -c darkblue --fastq ${lr}
    """
}


/*
 * STEP 1 - Read trimming and pre/post qc
 */
process fastqc_raw {
    tag "$name"
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") == -1 ? "QC_shortreads/fastqc/$filename" : null}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -t "${task.cpus}" -q $reads
    """
}


process fastp {
    tag "$name"
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "QC_shortreads/fastp/$name/$filename" : null}

    input:
    set val(name), file(reads) from read_files_fastp
    val adapter from params.adapter_forward
    val adapter_reverse from params.adapter_reverse
    val qual from params.mean_quality
    val trim_qual from params.trimming_quality

    output:
    set val(name), file("${name}_trimmed*.fastq.gz") into trimmed_reads
    file("fastp.*")

    script:
    def pe_input = params.singleEnd ? '' :  "-I \"${reads[1]}\""
    def pe_output1 = params.singleEnd ? "-o \"${name}_trimmed.fastq.gz\"" :  "-o \"${name}_trimmed_R1.fastq.gz\""
    def pe_output2 = params.singleEnd ? '' :  "-O \"${name}_trimmed_R2.fastq.gz\""
    """
    fastp -w "${task.cpus}" -q "${qual}" --cut_by_quality5 \
        --cut_by_quality3 --cut_mean_quality "${trim_qual}"\
        --adapter_sequence=${adapter} --adapter_sequence_r2=${adapter_reverse} \
        -i "${reads[0]}" $pe_input $pe_output1 $pe_output2
    """
}

/*
 * Remove PhiX contamination from Illumina reads
 * TODO: PhiX into/from iGenomes.conf?
 */
if(!params.keep_phix) {
    process phix_download_db {
        tag "${genome}"

        input:
        file(genome) from file_phix_db

        output:
        set file(genome), file("ref*") into phix_db

        script:
        """
        bowtie2-build --threads "${task.cpus}" "${genome}" ref
        """
    }

    process remove_phix {
        tag "$name"

        publishDir "${params.outdir}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "QC_shortreads/remove_phix/$filename" : null}

        input:
        set val(name), file(reads), file(genome), file(db) from trimmed_reads.combine(phix_db)

        output:
        set val(name), file("*.fastq.gz") into (trimmed_reads_megahit, trimmed_reads_metabat, trimmed_reads_fastqc, trimmed_sr_spadeshybrid, trimmed_reads_spades, trimmed_reads_centrifuge, trimmed_reads_kraken2, trimmed_reads_bowtie2)
        file("${name}_remove_phix_log.txt")

        script:
        if ( !params.singleEnd ) {
            """
            bowtie2 -p "${task.cpus}" -x ref -1 "${reads[0]}" -2 "${reads[1]}" --un-conc-gz ${name}_unmapped_%.fastq.gz
            echo "Bowtie2 reference: ${genome}" >${name}_remove_phix_log.txt
            zcat ${reads[0]} | echo "Read pairs before removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            zcat ${name}_unmapped_1.fastq.gz | echo "Read pairs after removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            """
        } else {
            """
            bowtie2 -p "${task.cpus}" -x ref -U ${reads}  --un-gz ${name}_unmapped.fastq.gz
            echo "Bowtie2 reference: ${genome}" >${name}_remove_phix_log.txt
            zcat ${reads[0]} | echo "Reads before removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            zcat ${name}_unmapped.fastq.gz | echo "Reads after removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            """
        }

    }
} else {
    trimmed_reads.into {trimmed_reads_megahit; trimmed_reads_metabat; trimmed_reads_fastqc; trimmed_sr_spadeshybrid; trimmed_reads_spades; trimmed_reads_centrifuge}
}


process fastqc_trimmed {
    tag "$name"
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") == -1 ? "QC_shortreads/fastqc/$filename" : null}

    input:
    set val(name), file(reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results_trimmed

    script:
    """
    fastqc -t "${task.cpus}" -q ${reads}
    """
}

/*
 * STEP - Taxonomic information
 */

process centrifuge_db_preparation {
    input:
    file(db) from file_centrifuge_db

    output:
    set val("${db.toString().replace(".tar.gz", "")}"), file("*.cf") into centrifuge_database

    script:
    """
    tar -xf "${db}"
    """
}

trimmed_reads_centrifuge
    .combine(centrifuge_database)
    .set { centrifuge_input }

process centrifuge {
    tag "${name}-${db_name}"
    publishDir "${params.outdir}/Taxonomy/centrifuge/${name}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".krona") == -1 ? filename : null}

    input:
    set val(name), file(reads), val(db_name), file(db) from centrifuge_input

    output:
    set val("centrifuge"), val(name), file("results.krona") into centrifuge_to_krona
    file("report.txt")
    file("kreport.txt")

    script:
    def input = params.singleEnd ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    centrifuge -x "${db_name}" \
        -p "${task.cpus}" \
        --report-file report.txt \
        -S results.txt \
        $input
    centrifuge-kreport -x "${db_name}" results.txt > kreport.txt
    cat results.txt | cut -f 1,3 > results.krona
    """
}

process kraken2_db_preparation {
    input:
    file(db) from file_kraken2_db

    output:
    set val("${db.baseName}"), file("${db.baseName}/*.k2d") into kraken2_database

    script:
    """
    tar -xf "${db}"
    """
}

trimmed_reads_kraken2
    .combine(kraken2_database)
    .set { kraken2_input }

process kraken2 {
    tag "${name}-${db_name}"
    publishDir "${params.outdir}/Taxonomy/kraken2/${name}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".krona") == -1 ? filename : null}

    input:
    set val(name), file(reads), val(db_name), file("database/*") from kraken2_input

    output:
    set val("kraken2"), val(name), file("results.krona") into kraken2_to_krona
    file("kraken2_report.txt")

    script:
    def input = params.singleEnd ? "\"${reads}\"" :  "--paired \"${reads[0]}\" \"${reads[1]}\""
    """
    kraken2 \
        --report-zero-counts \
        --threads "${task.cpus}" \
        --db database \
        --fastq-input \
        --report kraken2_report.txt \
        $input \
        > kraken2.kraken
    cat kraken2.kraken | cut -f 2,3 > results.krona
    """
}

process krona_db {
    output:
    file("taxonomy/taxonomy.tab") into file_krona_db

    when:
    ( params.centrifuge_db || params.kraken2_db ) && !params.skip_krona

    script:
    """
    ktUpdateTaxonomy.sh taxonomy
    """
}

centrifuge_to_krona
    .mix(kraken2_to_krona)
    .combine(file_krona_db)
    .set { krona_input }

process krona {
    tag "${classifier}-${name}"
    publishDir "${params.outdir}/Taxonomy/${classifier}/${name}", mode: 'copy'

    input:
    set val(classifier), val(name), file(report), file("taxonomy/taxonomy.tab") from krona_input

    output:
    file("*.html")

    script:
    """
    ktImportTaxonomy "$report" -tax taxonomy
    """
}


/*
 * STEP 2 - Assembly
 */
process megahit {
    tag "$name"
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "Assembly/$filename" : null}

    input:
    set val(name), file(reads) from trimmed_reads_megahit

    output:
    set val("MEGAHIT"), val("$name"), file("MEGAHIT/${name}.contigs.fa") into assembly_megahit_to_quast
    set val("MEGAHIT"), val("$name"), file("MEGAHIT/${name}.contigs.fa") into assembly_megahit_to_metabat
    file("MEGAHIT/*.log")

    when:
    !params.skip_megahit

    script:
    def input = params.singleEnd ? "-r \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    megahit -t "${task.cpus}" $input -o MEGAHIT --out-prefix "${name}"
    """
}


/*
 * metaSpades hybrid Assembly
 */

 files_lr_filtered
    .combine(trimmed_sr_spadeshybrid, by: 0)
    .set { files_pre_spadeshybrid }

process spadeshybrid {
    tag "$id"
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${id}*",
        saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "Assembly/SPAdesHybrid/$filename" : null}

    input:
    set id, file(lr), file(sr) from files_pre_spadeshybrid

    output:
    set id, val("SPAdesHybrid"), file("${id}_graph.gfa") into assembly_graph_spadeshybrid
    set val("SPAdesHybrid"), val("$id"), file("${id}_scaffolds.fasta") into assembly_spadeshybrid_to_quast
    set val("SPAdesHybrid"), val("$id"), file("${id}_scaffolds.fasta") into assembly_spadeshybrid_to_metabat
    file("${id}_contigs.fasta")
    file("${id}_log.txt")

    when:
    params.manifest && !params.singleEnd && !params.skip_spadeshybrid

    script:
    def maxmem = "${task.memory.toString().replaceAll(/[\sGB]/,'')}"
    """
    metaspades.py \
        --threads "${task.cpus}" \
        --memory "$maxmem" \
        --pe1-1 ${sr[0]} \
        --pe1-2 ${sr[1]} \
        --nanopore ${lr} \
        -o spades
    mv spades/assembly_graph_with_scaffolds.gfa ${id}_graph.gfa
    mv spades/scaffolds.fasta ${id}_scaffolds.fasta
    mv spades/contigs.fasta ${id}_contigs.fasta
    mv spades/spades.log ${id}_log.txt
    """
}


process spades {
    tag "$id"
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${id}*",
        saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "Assembly/SPAdes/$filename" : null}

    input:
    set id, file(sr) from trimmed_reads_spades

    output:
    set id, val("SPAdes"), file("${id}_graph.gfa") into assembly_graph_spades
    set val("SPAdes"), val("$id"), file("${id}_scaffolds.fasta") into assembly_spades_to_quast
    set val("SPAdes"), val("$id"), file("${id}_scaffolds.fasta") into assembly_spades_to_metabat
    file("${id}_contigs.fasta")
    file("${id}_log.txt")

    when:
    !params.singleEnd && !params.skip_spades

    script:
    def maxmem = "${task.memory.toString().replaceAll(/[\sGB]/,'')}"
    """
    metaspades.py \
        --threads "${task.cpus}" \
        --memory "$maxmem" \
        --pe1-1 ${sr[0]} \
        --pe1-2 ${sr[1]} \
        -o spades
    mv spades/assembly_graph_with_scaffolds.gfa ${id}_graph.gfa
    mv spades/scaffolds.fasta ${id}_scaffolds.fasta
    mv spades/contigs.fasta ${id}_contigs.fasta
    mv spades/spades.log ${id}_log.txt
    """
}


process quast {
    tag "$assembler-$sample"
    publishDir "${params.outdir}/Assembly/$assembler", mode: 'copy'

    input:
    set val(assembler), val(sample), file(assembly) from assembly_spades_to_quast.mix(assembly_megahit_to_quast).mix(assembly_spadeshybrid_to_quast)

    output:
    file("${sample}_QC/*") into quast_results

    when:
    !params.skip_quast

    script:
    """
    metaquast.py --threads "${task.cpus}" --rna-finding --max-ref-number 0 -l "${assembler}-${sample}" "${assembly}" -o "${sample}_QC"
    """
}

bowtie2_input = Channel.empty()

assembly_all_to_metabat = assembly_spades_to_metabat.mix(assembly_megahit_to_metabat,assembly_spadeshybrid_to_metabat)

(assembly_all_to_metabat, assembly_all_to_metabat_copy) = assembly_all_to_metabat.into(2)

bowtie2_input = assembly_all_to_metabat
    .combine(trimmed_reads_bowtie2)

(bowtie2_input, bowtie2_input_copy) = bowtie2_input.into(2)

/*
 * STEP 3 - Binning
 */
process bowtie2 {
    tag "$assembler-$sample"

    input:
    set val(assembler), val(sample), file(assembly), val(sampleToMap), file(reads) from bowtie2_input

    output:
    set val(assembler), val(sample), file("${assembler}-${sample}-${sampleToMap}.bam"), file("${assembler}-${sample}-${sampleToMap}.bam.bai") into assembly_mapping_for_metabat

    when:
    !params.skip_binning

    script:
    def name = "${assembler}-${sample}-${sampleToMap}"
    def input = params.singleEnd ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
        """
        bowtie2-build --threads "${task.cpus}" "${assembly}" ref
        bowtie2 -p "${task.cpus}" -x ref $input | \
            samtools view -@ "${task.cpus}" -bS | \
            samtools sort -@ "${task.cpus}" -o "${name}.bam"
        samtools index "${name}.bam"
        """
}

assembly_mapping_for_metabat = assembly_mapping_for_metabat.groupTuple(by:[0,1]).join(assembly_all_to_metabat_copy)

assembly_mapping_for_metabat = assembly_mapping_for_metabat.dump(tag:'assembly_mapping_for_metabat')

process metabat {
    tag "$assembler-$sample"
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> (filename.indexOf(".bam") == -1 && filename.indexOf(".fastq.gz") == -1) ? "GenomeBinning/$filename" : null}

    input:
    set val(assembler), val(sample), file(bam), file(index), val(sampleCopy), file(assembly) from assembly_mapping_for_metabat
    val(min_size) from params.min_contig_size

    output:
    set val(assembler), val(sample), file("MetaBAT2/*") into metabat_bins mode flatten
    set val(assembler), val(sample), file("MetaBAT2/*") into metabat_bins_for_cat
    set val(assembler), val(sample), file("MetaBAT2/*") into metabat_bins_quast_bins

    when:
    !params.skip_binning

    script:
    def name = "${assembler}-${sample}"
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o "MetaBAT2/${name}" -m ${min_size}

    #if bin folder is empty
    if [ -z \"\$(ls -A MetaBAT2)\" ]; then
        cp ${assembly} MetaBAT2/${assembler}-${assembly}
    fi
    """
}


process busco_download_db {
    tag "${database.baseName}"

    input:
    file(database) from file_busco_db

    output:
    set val("${database.toString().replace(".tar.gz", "")}"), file("buscodb/*") into busco_db

    script:
    """
    mkdir buscodb
    tar -xf ${database} -C buscodb
    """
}

metabat_bins
    .combine(busco_db)
    .set { metabat_db_busco }

/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */
process busco {
    tag "${assembly}"
    publishDir "${params.outdir}/GenomeBinning/QC/BUSCO/", mode: 'copy'

    input:
    set val(assembler), val(sample), file(assembly), val(db_name), file(db) from metabat_db_busco

    output:
    file("short_summary_${assembly}.txt") into (busco_summary_to_multiqc, busco_summary_to_plot)
    val("$assembler-$sample") into busco_assembler_sample_to_plot
    file("${assembly}_busco_log.txt")
    file("${assembly}_buscos.faa")
    file("${assembly}_buscos.fna")

    script:
    if( workflow.profile.toString().indexOf("conda") == -1) {
        """
        cp -r /opt/conda/pkgs/augustus*/config augustus_config/
        export AUGUSTUS_CONFIG_PATH=augustus_config

        run_BUSCO.py \
            --in ${assembly} \
            --lineage_path $db_name \
            --cpu "${task.cpus}" \
            --blast_single_core \
            --mode genome \
            --out ${assembly} \
            >${assembly}_busco_log.txt
        cp run_${assembly}/short_summary_${assembly}.txt short_summary_${assembly}.txt

        for f in run_${assembly}/single_copy_busco_sequences/*faa; do
            [ -e "\$f" ] && cat run_${assembly}/single_copy_busco_sequences/*faa >${assembly}_buscos.faa || touch ${assembly}_buscos.faa
            break
        done
        for f in run_${assembly}/single_copy_busco_sequences/*fna; do
            [ -e "\$f" ] && cat run_${assembly}/single_copy_busco_sequences/*fna >${assembly}_buscos.fna || touch ${assembly}_buscos.fna
            break
        done
        """
    } else {
        """
        run_BUSCO.py \
            --in ${assembly} \
            --lineage_path $db_name \
            --cpu "${task.cpus}" \
            --blast_single_core \
            --mode genome \
            --out ${assembly} \
            >${assembly}_busco_log.txt
        cp run_${assembly}/short_summary_${assembly}.txt short_summary_${assembly}.txt

        for f in run_${assembly}/single_copy_busco_sequences/*faa; do
            [ -e "\$f" ] && cat run_${assembly}/single_copy_busco_sequences/*faa >${assembly}_buscos.faa || touch ${assembly}_buscos.faa
            break
        done
        for f in run_${assembly}/single_copy_busco_sequences/*fna; do
            [ -e "\$f" ] && cat run_${assembly}/single_copy_busco_sequences/*fna >${assembly}_buscos.fna || touch ${assembly}_buscos.fna
            break
        done
        """
    }
}


process busco_plot {
    publishDir "${params.outdir}/GenomeBinning/QC/", mode: 'copy'

    input:
    file(summaries) from busco_summary_to_plot.collect()
    val(assemblersample) from busco_assembler_sample_to_plot.collect()

    output:
    file("*busco_figure.png")
    file("BUSCO/*busco_figure.R")
    file("BUSCO/*busco_summary.txt")
    file("busco_summary.txt") into busco_summary

    script:
    def assemblersampleunique = assemblersample.unique()
    """
    #for each assembler and sample:
    assemblersample=\$(echo \"$assemblersampleunique\" | sed 's/[][]//g')
    IFS=', ' read -r -a assemblersamples <<< \"\$assemblersample\"

    mkdir BUSCO

    for name in \"\${assemblersamples[@]}\"; do
        mkdir \${name}
        cp short_summary_\${name}* \${name}/
        generate_plot.py --working_directory \${name}

        cp \${name}/busco_figure.png \${name}-busco_figure.png
        cp \${name}/busco_figure.R \${name}-busco_figure.R

        summary_busco.py \${name}/short_summary_*.txt >BUSCO/\${name}-busco_summary.txt
    done

    cp *-busco_figure.R BUSCO/

    summary_busco.py short_summary_*.txt >busco_summary.txt
    """
}

process quast_bins {
    tag "$assembler-$sample"
    publishDir "${params.outdir}/GenomeBinning/QC/", mode: 'copy'

    input:
    set val(assembler), val(sample), file(assembly) from metabat_bins_quast_bins

    output:
    file("QUAST/*")
    file("QUAST/*-quast_summary.tsv") into quast_bin_summaries

    when:
    !params.skip_quast

    script:
    """
    ASSEMBLIES=\$(echo \"$assembly\" | sed 's/[][]//g')
    IFS=', ' read -r -a assemblies <<< \"\$ASSEMBLIES\"

    for assembly in \"\${assemblies[@]}\"; do
        metaquast.py --threads "${task.cpus}" --max-ref-number 0 --rna-finding --gene-finding -l "\${assembly}" "\${assembly}" -o "QUAST/\${assembly}"
        if ! [ -f "QUAST/${assembler}-${sample}-quast_summary.tsv" ]; then
            cp "QUAST/\${assembly}/transposed_report.tsv" "QUAST/${assembler}-${sample}-quast_summary.tsv"
        else
            tail -n +2 "QUAST/\${assembly}/transposed_report.tsv" >> "QUAST/${assembler}-${sample}-quast_summary.tsv"
        fi
    done
    """
}

process merge_quast_and_busco {
    publishDir "${params.outdir}/GenomeBinning/QC/", mode: 'copy'

    input:
    file(quast_bin_sum) from quast_bin_summaries.collect()
    file(busco_sum) from busco_summary

    output:
    file("quast_and_busco_summary.tsv")
    file("quast_summary.tsv")

    script:
    """
    QUAST_BIN=\$(echo \"$quast_bin_sum\" | sed 's/[][]//g')
    IFS=', ' read -r -a quast_bin <<< \"\$QUAST_BIN\"

    for quast_file in \"\${quast_bin[@]}\"; do
        if ! [ -f "quast_summary.tsv" ]; then
            cp "\${quast_file}" "quast_summary.tsv"
        else
            tail -n +2 "\${quast_file}" >> "quast_summary.tsv"
        fi
    done

    combine_tables.py $busco_sum quast_summary.tsv >quast_and_busco_summary.tsv
    """
}

/*
 * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
 */
process cat_db {
    tag "${database.baseName}"

    input:
    file(database) from file_cat_db

    output:
    set val("${database.toString().replace(".tar.gz", "")}"), file("database/*"), file("taxonomy/*") into cat_db

    script:
    """
    mkdir catDB
    tar -xf ${database} -C catDB
    mv `find catDB/ -type d -name "*taxonomy*"` taxonomy/
    mv `find catDB/ -type d -name "*CAT_database*"` database/
    """
}

metabat_bins_for_cat
    .combine(cat_db)
    .set { cat_input }

process cat {
    tag "${assembler}-${sample}-${db_name}"
    publishDir "${params.outdir}/Taxonomy/${assembler}", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".names.txt") > 0) filename
        else "raw/$filename"
    }

    input:
    set val(assembler), val(sample), file("bins/*"), val(db_name), file("database/*"), file("taxonomy/*") from cat_input

    output:
    file("*.ORF2LCA.txt")
    file("*.names.txt")
    file("*.predicted_proteins.faa")
    file("*.predicted_proteins.gff")
    file("*.log")
    file("*.bin2classification.txt")

    script:
    """
    CAT bins -b "bins/" -d database/ -t taxonomy/ -n "${task.cpus}" -s .fa --top 6 -o "${assembler}-${sample}" --I_know_what_Im_doing
    CAT add_names -i "${assembler}-${sample}.ORF2LCA.txt" -o "${assembler}-${sample}.ORF2LCA.names.txt" -t taxonomy/
    CAT add_names -i "${assembler}-${sample}.bin2classification.txt" -o "${assembler}-${sample}.bin2classification.names.txt" -t taxonomy/
    """
}

/*
 * STEP 4 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file (fastqc_raw:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file (fastqc_trimmed:'fastqc/*') from fastqc_results_trimmed.collect().ifEmpty([])
    file ('quast*/*') from quast_results.collect()
    //file ('short_summary_*.txt') from busco_summary_to_multiqc.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f ${rtitle} ${rfilename} --config ${multiqc_config} .
    """
}


/*
 * STEP 5 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/mag] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/mag] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/mag] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/mag] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/mag] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/mag] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/mag]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/mag]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/mag v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
