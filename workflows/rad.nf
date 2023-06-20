/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRad.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.genbank_ref) { ch_genbank_ref = file(params.genbank_ref) } else { exit 1, 'Genbank reference file not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MUGSY } from '../modules/local/mugsy'
include { MAKE_REFERENCE } from '../modules/local/make_reference'
include { GENERATE_CONSENSUS } from '../modules/local/generate_consensus'
include { GENBANK_TO_FASTA } from '../modules/local/genbank_to_fasta'
include { PROKKA_GENBANK_TO_FASTA_DB } from '../modules/local/prokka_genbank_to_fasta_db'
include { SUMMARY } from '../modules/local/summary'
include { CLEANUP } from '../modules/local/cleanup'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTQ_TRIM_FASTP_FASTQC     } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { FASTQ_ALIGN_BOWTIE2		  } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { FASTQ_ALIGN_BOWTIE2 as FASTQ_ALIGN_BOWTIE2_NEW_REF } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
// include { FASTQC                      } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SPADES } from '../modules/nf-core/spades/main' 
include { GUNZIP } from '../modules/nf-core/gunzip/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_ALIGNED } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED } from '../modules/nf-core/samtools/sort/main'
include { LAST_MAFCONVERT } from '../modules/nf-core/last/mafconvert/main'
include { GUNZIP as GUNZIP_MAF_TO_SAM } from '../modules/nf-core/gunzip/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_NEW_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { PROKKA } from '../modules/nf-core/prokka/main' 
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAD {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    GENBANK_TO_FASTA (
        ch_genbank_ref
    )
    
    PROKKA_GENBANK_TO_FASTA_DB (
        ch_genbank_ref
    )

    BOWTIE2_BUILD_REFERENCE (
        Channel.of([:]).combine(GENBANK_TO_FASTA.out.fasta)
    )

    ch_bowtie2_index = BOWTIE2_BUILD_REFERENCE.out.index.map{ [it[1]] }

    FASTQ_TRIM_FASTP_FASTQC (
        INPUT_CHECK.out.reads,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )
    
	FASTQ_ALIGN_BOWTIE2 ( 
        FASTQ_TRIM_FASTP_FASTQC.out.reads,
        FASTQ_TRIM_FASTP_FASTQC.out.reads.map { [it[0]] }.combine(ch_bowtie2_index),
		params.save_bowtie2_unaligned,
		params.sort_bowtie2_bam,
        FASTQ_TRIM_FASTP_FASTQC.out.reads.map { [it[0]] }.combine(GENBANK_TO_FASTA.out.fasta)
	)

    SPADES (
        FASTQ_TRIM_FASTP_FASTQC.out.reads.map { [it[0], it[1]]}
    )

    GUNZIP (
        SPADES.out.scaffolds
    )

    MUGSY (
        GUNZIP.out.gunzip.map {[it[0], it[1]]},
        GENBANK_TO_FASTA.out.fasta
    )

    LAST_MAFCONVERT (
        MUGSY.out.aligned_maf,
        "sam"
    )

    GUNZIP_MAF_TO_SAM (
        LAST_MAFCONVERT.out.sam_gz
    )
    
    SAMTOOLS_VIEW_ALIGNED (
        GUNZIP_MAF_TO_SAM.out.gunzip.map { [it[0], it[1], []] },
        GUNZIP_MAF_TO_SAM.out.gunzip.map { [it[0]] }.combine(GENBANK_TO_FASTA.out.fasta),
        []
    )

    SAMTOOLS_SORT_ALIGNED (
        SAMTOOLS_VIEW_ALIGNED.out.bam
    )

    MAKE_REFERENCE (
        SAMTOOLS_SORT_ALIGNED.out.bam,
        SAMTOOLS_SORT_ALIGNED.out.bam.map { [it[0]] }.combine(GENBANK_TO_FASTA.out.fasta),
        params.make_reference
    )

    BOWTIE2_BUILD_NEW_REFERENCE (
        MAKE_REFERENCE.out.new_ref
    )

	FASTQ_ALIGN_BOWTIE2_NEW_REF ( 
        FASTQ_TRIM_FASTP_FASTQC.out.reads,
        FASTQ_TRIM_FASTP_FASTQC.out.reads.map { [it[0]] }.join(BOWTIE2_BUILD_NEW_REFERENCE.out.index),
		params.save_bowtie2_unaligned,
		params.sort_bowtie2_bam,
        FASTQ_TRIM_FASTP_FASTQC.out.reads.map { [it[0]] }.join(MAKE_REFERENCE.out.new_ref)
	)

    GENERATE_CONSENSUS (
        FASTQ_ALIGN_BOWTIE2_NEW_REF.out.bam,
        FASTQ_ALIGN_BOWTIE2_NEW_REF.out.bai,
        FASTQ_ALIGN_BOWTIE2_NEW_REF.out.bam.map { [it[0]] }.combine(GENBANK_TO_FASTA.out.fasta),
        params.generate_consensus
    )

    PROKKA (
        GENERATE_CONSENSUS.out.consensus,
        PROKKA_GENBANK_TO_FASTA_DB.out.faa,
        []
    )

    FASTQ_TRIM_FASTP_FASTQC.out.trim_log
        .join(FASTQ_ALIGN_BOWTIE2_NEW_REF.out.bam)
        .join(FASTQ_ALIGN_BOWTIE2_NEW_REF.out.bai)
        .join(GENERATE_CONSENSUS.out.consensus)
        .map { meta, trim_log, bam, bai, consensus -> [ meta, trim_log, bam, bai, consensus] }
        .set {ch_summary}
    
    SUMMARY (
        ch_summary
    )
        
    CLEANUP ( SUMMARY.out.summary_tsv.collect() )    

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowCmv.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowCmv.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/