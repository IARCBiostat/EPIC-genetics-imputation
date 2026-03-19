#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { NORMALIZATION } from './modules/normalization.nf'
include { LIFTOVER } from './modules/liftover.nf'
include { PHASING } from './modules/phasing.nf'
include { IMPUTATION } from './modules/imputation.nf'
include { POST_IMPUTATION_VARIANT_QC; POST_IMPUTATION_SAMPLE_QC } from './modules/post_imputation.nf'

workflow STAGE_NORMALIZATION {
    take: ch_studies
    main:
        NORMALIZATION(ch_studies)
    emit:
        normalized = NORMALIZATION.out.normalized_data
}

workflow STAGE_LIFTOVER {
    take: ch_normalized
    main:
        ch_chain = Channel.fromPath(params.chain_file)
        LIFTOVER(ch_normalized.combine(ch_chain))
    emit:
        lifted = LIFTOVER.out.lifted_data
}

workflow STAGE_IMPUTATION {
    take: ch_lifted
    main:
        ch_chrs = Channel.fromList(params.chromosomes)
        ch_vcf_prep = ch_lifted
            .combine(ch_chrs)
            .map { sid, bed, bim, fam, chr -> 
                tuple(sid, chr, file("${sid}_chr${chr}.vcf.gz")) 
            }

        ch_maps = Channel.fromPath("${params.genetic_maps_dir}/chr*.b37.gmap.gz")
            .map { f -> tuple(f.name.replaceAll(/chr(.*)\.b37\.gmap\.gz/, '$1'), f) }
        
        ch_phasing_in = ch_vcf_prep.join(ch_maps, by: 1)
            .map { chr, sid, vcf, map -> tuple(sid, chr, vcf, map) }

        PHASING(ch_phasing_in)

        ch_ref_vcfs = Channel.fromPath("${params.ref_1000g_dir}/ALL.chr*.phase3*.vcf.gz")
            .map { f -> tuple(f.name.replaceAll(/ALL\.chr(.*)\.phase3.*/, '$1'), f) }
        
        ch_impute_in = PHASING.out.phased_vcf.join(ch_ref_vcfs, by: 1)
            .map { chr, sid, phased, ref -> tuple(sid, chr, phased, ref) }

        IMPUTATION(ch_impute_in)
    emit:
        imputed = IMPUTATION.out.imputed_vcf
}

workflow {
    log.info """
        EPIC IMPUTATION PIPELINE
        ========================
        samplesheet : ${params.samplesheet}
        outdir      : ${params.outdir}
        """

    // 1. Validation
    VALIDATE_INPUTS(file(params.samplesheet))
    
    // 2. Load samplesheet
    ch_studies = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .filter { row -> params.include_studies == "" || params.include_studies.split(',').contains(row.study_id) }
        .map { row -> 
            tuple(row.study_id, file(row.bed), file(row.bim), file(row.fam), file(row.manifest), file(row.id_linkage))
        }

    // Run the full pipeline
    STAGE_NORMALIZATION(ch_studies)
    STAGE_LIFTOVER(STAGE_NORMALIZATION.out.normalized)
    STAGE_IMPUTATION(STAGE_LIFTOVER.out.lifted)
    
    POST_IMPUTATION_VARIANT_QC(STAGE_IMPUTATION.out.imputed)
    
    ch_sample_qc_in = POST_IMPUTATION_VARIANT_QC.out.filtered_plink
        .groupTuple(by: 0)
        .map { sid, chrs, beds, bims, fams -> tuple(sid, beds, bims, fams) }
    
    POST_IMPUTATION_SAMPLE_QC(ch_sample_qc_in)
}
