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
if (params.bowtie2_host_index) { ch_bowtie2_host_index = Channel.fromPath(params.bowtie2_host_index)} else { ch_bowtie2_host_index = [] }
if (params.region_map) { ch_region_map = file(params.region_map) } else { ch_region_map = [] }

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
include { SRA_SCRUB_HUMAN } from '../modules/local/sra_scrub_human.nf'
include { BWA_MEM_ALIGN as BWA_MEM_ALIGN_NEW_REF } from '../modules/local/bwa_mem_align'
include { IVAR_CONSENSUS } from '../modules/local/ivar_consensus'
include { IVAR_VARIANTS } from '../modules/local/ivar_variants'
include { FORMAT_VARIANTS } from '../modules/local/format_variants'
include { BWA_MEM_ALIGN as BWA_MEM_ALIGN_FINAL_CONSENSUS } from '../modules/local/bwa_mem_align'
include { SUMMARY } from '../modules/local/summary'
include { CLEANUP } from '../modules/local/cleanup'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
//include { FASTQ_TRIM_FASTP_FASTQC     } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
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
include { UNICYCLER } from '../modules/nf-core/unicycler/main'
include { GUNZIP } from '../modules/nf-core/gunzip/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_ALIGNED } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED } from '../modules/nf-core/samtools/sort/main'
include { LAST_MAFCONVERT } from '../modules/nf-core/last/mafconvert/main'
include { GUNZIP as GUNZIP_MAF_TO_SAM } from '../modules/nf-core/gunzip/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_NEW_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_FINAL_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_HOST_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_NEW_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_FINAL_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { BBMAP_BBDUK as BBDUK_R } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_L } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_Q } from '../modules/nf-core/bbmap/bbduk/main'
include { FASTQC as FASTQC_RAW } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/fastqc/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAD {

    ch_versions = Channel.empty()

    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    GENBANK_TO_FASTA (
        ch_genbank_ref,
        ch_region_map,
        params.configure_reference
    )

    BOWTIE2_BUILD_REFERENCE (
        Channel.of([:]).combine(GENBANK_TO_FASTA.out.fasta)
    )

    ch_bowtie2_index = BOWTIE2_BUILD_REFERENCE.out.index.map{ [it[1]] }

    
    ch_raw_reads = INPUT_CHECK.out.reads

    if (params.sub_sample) {
        SEQTK_SAMPLE (
            ch_raw_reads.combine([params.sample_size])
        )
        ch_raw_reads = SEQTK_SAMPLE.out.reads
    }

    FASTQC_RAW (
        ch_raw_reads
    )

    BBDUK_R (
        ch_raw_reads,
        []
    )
    
    BBDUK_L (
        BBDUK_R.out.reads,
        []
    )

    BBDUK_Q (
        BBDUK_L.out.reads,
        []
    )

    ch_trimmed_reads = BBDUK_Q.out.reads

    if (params.bowtie2_host_index) {
    
        BOWTIE2_ALIGN_HOST_REFERENCE (
            BBDUK_Q.out.reads,
            BBDUK_Q.out.reads.map { [it[0]] }.combine(ch_bowtie2_host_index),
            true,  // save_unaligned, i.e. non-host reads
            false  // don't sort host aligned reads      
        )

        ch_trimmed_reads = BOWTIE2_ALIGN_HOST_REFERENCE.out.fastq
    }

    if (params.scrub_human) {
        SRA_SCRUB_HUMAN (
            ch_trimmed_reads
        )
        ch_trimmed_reads = SRA_SCRUB_HUMAN.out.reads
    }

	BOWTIE2_ALIGN_REFERENCE (
        ch_trimmed_reads,
        ch_trimmed_reads.map { [it[0]] }.combine(ch_bowtie2_index),
		params.save_bowtie2_unaligned,
		params.sort_bowtie2_bam
	)

    FASTQC_TRIMMED (
       ch_trimmed_reads
    )

    ch_scaffolds = Channel.empty()

    if (params.spades_flag == "unicycler") {
        
        UNICYCLER (
            ch_trimmed_reads
        )
        ch_scaffolds = UNICYCLER.out.scaffolds

    } else {
        
        SPADES (
            ch_trimmed_reads,
            params.spades_flag
        )
        ch_scaffolds = SPADES.out.scaffolds
    }

    GUNZIP (
        ch_scaffolds
    )

    if (params.region_map) { ch_regions = GENBANK_TO_FASTA.out.regions } else { ch_regions = [] }

    MUGSY (
        GUNZIP.out.gunzip.map {[it[0], it[1]]},
        GENBANK_TO_FASTA.out.fasta,
        ch_genbank_ref,
        ch_regions,
        params.maf_convert,
        params.make_reference,
        ch_region_map,
        params.configure_reference
    )

    BOWTIE2_BUILD_NEW_REFERENCE (
       MUGSY.out.new_ref
    )

    ch_trimmed_reads
        .join(BOWTIE2_BUILD_NEW_REFERENCE.out.index).set{ch_matched_reference}

    BOWTIE2_ALIGN_NEW_REFERENCE (
        ch_matched_reference.map{ [it[0], it[1]] },
        ch_matched_reference.map{ [it[0], it[2]] },
        false,
        true
    )

    IVAR_CONSENSUS (
        BOWTIE2_ALIGN_NEW_REFERENCE.out.bam
    )

    BOWTIE2_BUILD_FINAL_REFERENCE (
       IVAR_CONSENSUS.out.consensus
    )

    IVAR_VARIANTS (
        BOWTIE2_ALIGN_REFERENCE.out.bam,
        GENBANK_TO_FASTA.out.fasta,
        GENBANK_TO_FASTA.out.genes_gff
    )

    FORMAT_VARIANTS (
        IVAR_VARIANTS.out.variants,
        GENBANK_TO_FASTA.out.genes_gff,
        params.edit_ivar_variants
    )

    ch_trimmed_reads
        .join(BOWTIE2_BUILD_FINAL_REFERENCE.out.index).set{ch_reads_mapped_final_consensus}

    BOWTIE2_ALIGN_FINAL_REFERENCE (
        ch_reads_mapped_final_consensus.map{ [it[0], it[1]] },
        ch_reads_mapped_final_consensus.map{ [it[0], it[2]] },
        false,
        true
    )        

    BBDUK_R.out.log
        .join(BBDUK_Q.out.log)
        .join(BOWTIE2_ALIGN_NEW_REFERENCE.out.bam)
        .join(IVAR_CONSENSUS.out.consensus)
        .map { meta, rlog, qlog, bam, consensus -> [ meta, rlog, qlog, bam, consensus] }
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
