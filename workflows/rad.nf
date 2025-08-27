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
if (!params.kraken_host_db) { exit 1, "Kraken host database not specified with e.g. '--kraken_host_db s3://fh-pi-jerome-k-eco/greninger-lab/greninger-lab-file-share/refs/Kraken2_human/k2_human/' or via a detectable config file." }
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
include { BBMAP_REPAIR } from '../modules/local/bbmap_repair'
include { GENBANK_TO_FASTA } from '../modules/local/genbank_to_fasta'
include { FASTA_INDEX } from '../modules/local/fasta_index'
include { BBMAP_DEDUPE } from '../modules/local/bbmap_dedupe.nf'
include { READS_MAPPED } from '../modules/local/reads_mapped.nf'
include { PICARD_ADDORREPLACEREADGROUPS } from '../modules/local/addorreplacereadgroups'
include { GATK_REALIGNERTARGETCREATOR   } from '../modules/local/realignertargetcreator'
include { GATK_INDELREALIGNER           } from '../modules/local/indelrealigner'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_NEW_REF } from '../modules/local/addorreplacereadgroups'
include { GATK_REALIGNERTARGETCREATOR as GATK_REALIGNERTARGETCREATOR_NEW_REF   } from '../modules/local/realignertargetcreator'
include { GATK_INDELREALIGNER as GATK_INDELREALIGNER_NEW_REF           } from '../modules/local/indelrealigner'
include { GET_BEST_REFERENCE } from '../modules/local/get_best_reference'
include { IVAR_CONSENSUS } from '../modules/local/ivar_consensus'
include { IVAR_VARIANTS } from '../modules/local/ivar_variants'
include { FORMAT_VARIANTS } from '../modules/local/format_variants'
include { SUMMARY } from '../modules/local/summary'
include { SUMMARY as SUMMARY_FAILED } from '../modules/local/summary'
include { CLEANUP } from '../modules/local/cleanup'
include { CLEANUP as CLEANUP_FAILED } from '../modules/local/cleanup'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
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
include { SPADES } from '../modules/nf-core/spades/main'
include { UNICYCLER } from '../modules/nf-core/unicycler/main'
include { GUNZIP } from '../modules/nf-core/gunzip/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_NEW_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_FINAL_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_NEW_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_FINAL_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BBMAP_BBDUK as BBDUK_R } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_L } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_Q } from '../modules/nf-core/bbmap/bbduk/main'
include { FASTQC as FASTQC_RAW } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/fastqc/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAD {

    //ch_versions = Channel.empty()

    INPUT_CHECK (
        ch_input
    )

    GENBANK_TO_FASTA (
        ch_genbank_ref,
        ch_region_map,
        params.configure_reference
    )

    BOWTIE2_BUILD_REFERENCE (
        Channel.of([:]).combine(GENBANK_TO_FASTA.out.fasta)
    )

    ch_bowtie2_index = BOWTIE2_BUILD_REFERENCE.out.index.map{ [it[1]] }

    ch_raw_reads = INPUT_CHECK.out.reads.map{ [it[0], it[1]] }

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

    KRAKEN2_KRAKEN2 (
        BBDUK_Q.out.reads,
        params.kraken_host_db,
        true,  // save fastqs
        false  // don't report      
    )

    BBMAP_REPAIR (
        KRAKEN2_KRAKEN2.out.unclassified_reads_fastq,
    )    

    ch_trimmed_reads = BBMAP_REPAIR.out.fastqs

	BOWTIE2_ALIGN_REFERENCE (
        ch_trimmed_reads,
        ch_trimmed_reads.map { [it[0]] }.combine(ch_bowtie2_index),
		params.save_bowtie2_unaligned,
		params.sort_bowtie2_bam
	)

    FASTQC_TRIMMED (
        ch_trimmed_reads
    )

    PICARD_ADDORREPLACEREADGROUPS (
        BOWTIE2_ALIGN_REFERENCE.out.bam,
        [[],[]],
        [[],[]]
    )

    GATK_REALIGNERTARGETCREATOR (
        PICARD_ADDORREPLACEREADGROUPS.out.bam,
        PICARD_ADDORREPLACEREADGROUPS.out.bam.map{ [it[0]] }.combine(GENBANK_TO_FASTA.out.fasta),
        PICARD_ADDORREPLACEREADGROUPS.out.bam.map{ [it[0]] }.combine(GENBANK_TO_FASTA.out.fai),
        PICARD_ADDORREPLACEREADGROUPS.out.bam.map{ [it[0]] }.combine(GENBANK_TO_FASTA.out.dict),
        [[],[]]
    )

    GATK_INDELREALIGNER (
        PICARD_ADDORREPLACEREADGROUPS.out.bam.join(GATK_REALIGNERTARGETCREATOR.out.intervals),
        PICARD_ADDORREPLACEREADGROUPS.out.bam.map{ [it[0]] }.combine(GENBANK_TO_FASTA.out.fasta),
        PICARD_ADDORREPLACEREADGROUPS.out.bam.map{ [it[0]] }.combine(GENBANK_TO_FASTA.out.fai),
        PICARD_ADDORREPLACEREADGROUPS.out.bam.map{ [it[0]] }.combine(GENBANK_TO_FASTA.out.dict),
        [[],[]]
    )

    BBMAP_DEDUPE (
        ch_trimmed_reads
    )

    READS_MAPPED (
        BOWTIE2_ALIGN_REFERENCE.out.bam        
    )

    READS_MAPPED.out.reads_mapped.branch {
        passed: it[1].toInteger() > 50000
        failed: it[1].toInteger() <= 50000
    }.set{ ch_reads_mapped }

    BBMAP_DEDUPE.out.reads.join(ch_reads_mapped.passed).map{ [it[0], it[1]] }.set{ ch_trimmed_reads_passed }
    BBMAP_DEDUPE.out.reads.join(ch_reads_mapped.failed).map{ [it[0], it[1]] }.set{ ch_trimmed_reads_failed }

    ch_scaffolds = Channel.empty()

    if (params.spades_flag == "unicycler") {
        
        UNICYCLER (
            ch_trimmed_reads_passed
        )
        ch_scaffolds = UNICYCLER.out.scaffolds

    } else {
        
        SPADES (
            ch_trimmed_reads_passed,
            params.spades_flag
        )
        ch_scaffolds = SPADES.out.scaffolds
    }

    GUNZIP (
        ch_scaffolds
    )

    if (params.find_reference) {
        
        GUNZIP.out.gunzip.map {[it[0], it[1]]}
            .join(ch_trimmed_reads_passed)
            .join(BOWTIE2_ALIGN_REFERENCE.out.bam).set{ch_scaffolds}

        GET_BEST_REFERENCE (
            ch_scaffolds,
            params.find_reference
        )
            
        ch_scaffolds
            .join(GET_BEST_REFERENCE.out.regions)
            .join(GET_BEST_REFERENCE.out.region_map)
            .join(GET_BEST_REFERENCE.out.gb)
            //.map { meta, scaffolds, regions, region_map, gb -> [ meta, scaffolds, regions, region_map, gb ] }
            .set {ch_best_ref}

        MUGSY (
            ch_best_ref.map {[it[0], it[1]]},
            [],
            ch_best_ref.map {it[6]},
            ch_best_ref.map {it[4]},
            params.maf_convert,
            params.make_reference,
            ch_best_ref.map {[it[5]]},
            params.configure_reference
        )

    } else {

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
    }

    FASTA_INDEX (
        MUGSY.out.new_ref
    )    

    BOWTIE2_BUILD_NEW_REFERENCE (
        MUGSY.out.new_ref
    )

    ch_trimmed_reads_passed
        .join(BOWTIE2_BUILD_NEW_REFERENCE.out.index).set{ch_matched_reference}

    BOWTIE2_ALIGN_NEW_REFERENCE (
        ch_matched_reference.map{ [it[0], it[1]] },
        ch_matched_reference.map{ [it[0], it[2]] },
        false,
        true
    )

    PICARD_ADDORREPLACEREADGROUPS_NEW_REF (
        BOWTIE2_ALIGN_NEW_REFERENCE.out.bam,
        [[],[]],
        [[],[]]
    )

    PICARD_ADDORREPLACEREADGROUPS_NEW_REF.out.bam
    .join(FASTA_INDEX.out.fasta)
    .join(FASTA_INDEX.out.fai)
    .join(FASTA_INDEX.out.dict).set{ch_new_ref}

    GATK_REALIGNERTARGETCREATOR_NEW_REF (
        ch_new_ref.map{[it[0],it[1],it[2]]},
        ch_new_ref.map{[it[0],it[3]]},
        ch_new_ref.map{[it[0],it[4]]},
        ch_new_ref.map{[it[0],it[5]]},
        [[],[]]
    )

    PICARD_ADDORREPLACEREADGROUPS_NEW_REF.out.bam
    .join(GATK_REALIGNERTARGETCREATOR_NEW_REF.out.intervals)
    .join(FASTA_INDEX.out.fasta)
    .join(FASTA_INDEX.out.fai)
    .join(FASTA_INDEX.out.dict).set{ch_new_ref_interval}

    GATK_INDELREALIGNER_NEW_REF (
        ch_new_ref_interval.map{[it[0],it[1],it[2],it[3]]},
        ch_new_ref_interval.map{[it[0],it[4]]},
        ch_new_ref_interval.map{[it[0],it[5]]},
        ch_new_ref_interval.map{[it[0],it[6]]},
        [[],[]]
    )  

    IVAR_CONSENSUS (
        GATK_INDELREALIGNER_NEW_REF.out.bam.map{ [it[0], it[1]] }
    )

    BOWTIE2_BUILD_FINAL_REFERENCE (
       IVAR_CONSENSUS.out.consensus
    )

    IVAR_VARIANTS (
        GATK_INDELREALIGNER.out.bam.map{ [it[0], it[1]] },
        GENBANK_TO_FASTA.out.fasta,
        GENBANK_TO_FASTA.out.genes_gff
    )

    FORMAT_VARIANTS (
        IVAR_VARIANTS.out.variants,
        GENBANK_TO_FASTA.out.fasta,
        GENBANK_TO_FASTA.out.genes_gff,
        params.edit_ivar_variants
    )

    ch_trimmed_reads_passed
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
        .join(ch_trimmed_reads_passed)
        .map { meta, rlog, qlog, bam, consensus, reads_mapped -> [ meta, rlog, qlog, bam, consensus, reads_mapped] }
        .set {ch_summary_passed}

    BBDUK_R.out.log
        .join(BBDUK_Q.out.log)
        .join(BOWTIE2_ALIGN_REFERENCE.out.bam)
        .join(READS_MAPPED.out.empty)
        .join(ch_trimmed_reads_failed)
        .map { meta, rlog, qlog, bam, empty, reads_mapped -> [ meta, rlog, qlog, bam, empty, reads_mapped] }
        .set {ch_summary_failed}

    SUMMARY (
        ch_summary_passed
    )

    SUMMARY_FAILED (
        ch_summary_failed
    )    

    CLEANUP ( 
        SUMMARY.out.summary_tsv.collect(),
        false
    )    

    CLEANUP_FAILED (
        SUMMARY_FAILED.out.summary_tsv.collect(),
        true
    )

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
