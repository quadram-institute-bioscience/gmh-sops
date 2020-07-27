#!/usr/bin/env nextflow
//nextflow.preview.dsl=2
version = '0.1.2'

/*
 HISTORY
 0.1.2 -- Abricate summary MQC
*/

def helpMessage() {

    log.info"""
    Usage:
    nextflow run metagmh --reads '*_R{1,2}.fastq.gz' -profile docker

    Main arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --reference                   Reference genome
      --outdir                      Output directory 
      --compare                     Enable multisample comparisons (default: OFF)

    Optional arguments:  
      
      --tempdir                     Absolute PATH to the temporary directory (default ./tmp)
      --assembler                   Assembler to use. Default 'shovill' (optional: 'unicycler')
      --kraken2db                   Path to Kraken2 database
      --genus                       Genus (can be autodetected, for labeling only)
      --species                     Species (can be autodetected, for labeling only)

      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, singularity, test and more.

    """.stripIndent()
}
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

params.readPaths = false
params.compare = false
params.genus = '(genus)'
params.species = '(sp.)'
params.reads = ""

params.cpus = 6
params.assembler = 'shovill'
params.outdir = "nullarflow"
params.tempdir = "/tmp/"
params.adapter_forward  = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse  = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality     = 18
params.trimming_quality = 14

params.subsample        = 0.6
params.kraken2db        = '/qi/db/kraken2/mini/'

params.multiqc_config = "$baseDir/assets/multiqc.yaml"
Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }

if (params.readPaths) {
    if (params.singleEnd) {

        read_pairs_ch = Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            
    } else {

        read_pairs_ch = Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            
    }
} else {

    read_pairs_ch = Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        
}




log.info """
===        GMH - Nullarflow v.${version}          ===
============================================
 reads        : ${params.reads}
 tempdir      : ${params.tempdir}
 outdir       : ${params.outdir}
 """

 process versions {
   publishDir params.outdir, mode: "copy"
   label 'onecore'

   output:
   file('versions.txt') into versions_ch
   file('software_versions_mqc.html') into versions_multiqc_ch

   script:
   """
   software.pl --html software_versions_mqc.html > versions.txt
   """
 }


process filter_reads {
    publishDir "${params.outdir}/reads", mode: "copy"
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
    set val(name), file("${name}.filtered*.fastq.gz") into cleaned_unicycler_ch, cleaned_kraken_ch, cleaned_skesa_ch
    file("${name}.fastp.json") into json_ch

    script:
    if ( params.singleEnd )
        """
        fastp -w "${task.cpus}" -q "${qual}" --cut_by_quality5 \
          --cut_by_quality3 --cut_mean_quality "${trim_qual}" \
            --json ${name}.fastp.json \
            --adapter_sequence=${adapter} --adapter_sequence_r2=${adapter_reverse} \
            -i "${reads[0]}"  \
            -o "${name}.filtered_R1.fastq.gz"
        """
    else 
        """
        fastp -w "${task.cpus}" -q "${qual}" --cut_by_quality5 \
          --cut_by_quality3 --cut_mean_quality "${trim_qual}" \
            --json ${name}.fastp.json \
            --adapter_sequence=${adapter} --adapter_sequence_r2=${adapter_reverse} \
            -i "${reads[0]}" -I "${reads[1]}" \
            -o "${name}.filtered_R1.fastq.gz" -O "${name}.filtered_R2.fastq.gz"
        """

}

process assembly  {
    publishDir params.outdir, mode: "copy"
    tag "$name"

    input:
    set val(name), file(reads)   from cleaned_unicycler_ch

    output:
    set val(name), file("${name}.fa") into assembly_mlst_ch, assembly_abricate_ch,assembly_prokka_ch, quast_input_ch

    script:
    if ( params.assembler == 'shovill' )
        """
        shovill --assembler skesa --R1 "${reads[0]}" --R2 "${reads[1]}" --outdir assembly --cpus ${task.cpus} --ram 12.0 
        mv assembly/contigs.fa ${name}.fa
        """
    else if ( params.assembler == 'unicycler' )
        """
        unicycler -1 "${reads[0]}" -2 "${reads[1]}" --out assembly
        mv assembly/assembly.fasta ${name}.fa
        """
    else
        throw new IllegalArgumentException("Unknown assembler: $params.assembler")


}

process kraken2 {
    publishDir params.outdir, mode: "copy"
    tag "$name"

    input:
    set val(name), file(reads)   from cleaned_kraken_ch

    output:
    file("*.taxonomy.txt")       into kraken2_ch
    file("${name}.taxonomy.txt") into taxonomy_ch

    script:
    """
    kraken2 --threads ${task.cpus} --db ${params.kraken2db}  \
      --report ${name}.taxonomy.txt --paired --confidence 0.1 "${reads[0]}" "${reads[1]}" > /dev/null    
    """
}

process mlst {
    publishDir params.outdir, mode: "copy"
    tag "$name"

    input:
    set val(name), file('assembly.fasta')   from assembly_mlst_ch

    output:
    set val(name), file("*.txt") into mlst_ch

    script:
    """
    mlst assembly.fasta > ${name}.mlst.txt   
    """
}

process abricate {
    publishDir params.outdir, mode: "copy"
    tag "$name"

    input:
    set val(name), file("${name}")   from assembly_abricate_ch

    output:
    set val(name), file("${name}.*.tab") into abricate_ch 
    file("${name}*.tab") into abricate_summary_ch

    script:
    """
    abricate  --db resfinder --quiet --threads ${task.cpus} ${name} > ${name}.resistome.tab  
    abricate  --db vfdb      --quiet --threads ${task.cpus} ${name} > ${name}.virulome.tab   
    """
}

process abricate_summary {
    tag "report"
    publishDir "${params.outdir}", mode: "copy"

    input:
    file(files) from abricate_summary_ch.collect().ifEmpty([])

    output:
    file('virulome_mqc.txt')  into virulome_mqc
    file('resistome_mqc.txt') into resistome_mqc

    script:
    """
    abricate --summary *.virulome.tab  | sed 's/.virulome.tab//'  > virulome_mqc.txt 
    abricate --summary *.resistome.tab | sed 's/.resistome.tab//' > resistome_mqc.txt 
    """
}

process prokka {
    publishDir params.outdir, mode: "copy"
    tag "$name"

    input:
    set val(name), file('contigs.fa')   from assembly_prokka_ch
    file('taxon.txt')  from taxonomy_ch

    output:
    set val(name), file("prokka/*.faa") into proteins_ch
    file("prokka/*.gff") into annotation_ch
    file("${name}.txt") into prokka_ch

    script:
    """
    # 
    prokka `class_to_param.pl taxon.txt` --outdir prokka  --locustag $name --strain $name --prefix $name --cpus ${task.cpus} contigs.fa
    mv prokka/*.txt ${name}.txt
    """
}

process roary {
    tag 'report'
    publishDir params.outdir, mode: "copy"

    input:
    file(annotations) from annotation_ch.collect().ifEmpty([])

    output:
    file('roary.tsv')
    file('roary.zip')

    script:
    """
    mkdir roary
    cd roary
    roary -e --mafft -p ${task.cpus}  ../*.gff -o roary
    mv roary ../roary.tsv
    
    if [ \$(command -v tar) ]; then
        tar cvfz ../roary.tar.gz *.*
    fi
    """
}

process quast {
    publishDir "${params.outdir}/", mode: "copy"

    input:
    file('*')   from quast_input_ch.collect().ifEmpty([])

    output:
    file('quast') into quast_dir_ch
    file('quast/report.tsv') into quast_qc_ch

    script:
    """
    quast -o quast *.fa
    """
}

process multiqc {
    publishDir params.outdir, mode: "copy"

    input:
    file("fastp/*")         from json_ch.collect().ifEmpty([])
    file("kraken2/*")       from kraken2_ch.collect().ifEmpty([])
    file("*.txt")           from prokka_ch.collect().ifEmpty([])
    file multiqc_config     from ch_config_for_multiqc
    file('virulome_mqc.tsv')  from virulome_mqc
    file('resistome_mqc.tsv') from resistome_mqc
    file('report.tsv')      from quast_qc_ch
    file('software_versions_mqc.html') from versions_multiqc_ch

    output:
    file "*.html" into multiqc_report
    file "*_data"


    script:
    reportTitle    = custom_runName ? "--title \"$custom_runName\"" : ''
    reportFilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    sed -i 's/_R1//g' fastp/*.json
    for i in kraken2/*.txt; do
      mv \$i \${i/.taxonomy.txt/};
    done
    multiqc --config $multiqc_config --comment "Automatic output from Nullarflow" $reportTitle $reportFilename . 
    """
}
workflow.onComplete {
    log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    log.info ( workflow.success ?
        "\nDone! The results are saved in --> $params.outdir/\n" :
        "Oops .. something went wrong" )
}
