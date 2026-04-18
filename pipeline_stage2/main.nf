#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * MAIN EPIC IMPUTATION PIPELINE
 * Input datasets come from stage 1 outputs in analysis/<STUDY>/stage1/
 */

include { PREP_REFERENCE } from './modules/reference_prep.nf'
include { PREP_TARGET_VCF } from './modules/target_prep.nf'
include { PHASE_AUTOSOMES } from './modules/phasing.nf'
include { IMPUTE_AUTOSOMES } from './modules/imputation.nf'
include { IMPUTE_CHRX } from './modules/chrX_imputation.nf'

workflow {

    // Define standard chromosomes
    ch_chroms_auto = Channel.fromList((1..22).collect { it.toString() })
    ch_chroms_X = Channel.of('X')

    // Parse included studies
    def included_studies = params.study == 'all' ? [] : params.study.split(',').collect{ it.trim() }

    // Read final hg38 PLINK files from analysis/<study>/stage1/<study>.{bed,bim,fam}
    ch_stage1_plink = Channel.fromFilePairs("${params.stage1_root}/*/stage1/*.{bed,bim,fam}", size: 3, flat: true)
        .map { prefix, bed, bim, fam -> tuple(bed.baseName, prefix, bed, bim, fam) }
        .filter { study_name, prefix, bed, bim, fam -> 
            params.study == 'all' || included_studies.contains(study_name) 
        }

    // 1. Reference pre-processing and stage-1 target export
    if (params.run_preprocessing) {
        ch_reference = PREP_REFERENCE(ch_chroms_auto.mix(ch_chroms_X))

        // Convert stage-1 hg38 PLINK inputs to one indexed target VCF per chromosome
        ch_target_vcf = PREP_TARGET_VCF(ch_stage1_plink)

        // Flatten chromosome-specific target VCF and index pairs
        ch_target_vcf_flat = ch_target_vcf.transpose()
            .map { study_name, vcf, tbi ->
                def matcher = vcf.name =~ /_chr([0-9X]+)\.vcf\.gz$/
                def chr = matcher[0][1]
                tuple(chr, study_name, vcf, tbi)
            }
            // Use combine-by-chromosome rather than join so every study/chromosome
            // target VCF receives the shared reference for that chromosome.
            .combine(ch_reference, by: 0)
            .branch {
                autosomes: it[0] != 'X'
                chrX: it[0] == 'X'
            }
    }

    // 2. Phasing using Eagle v2.4.1
    if (params.run_phasing) {
        // [chr, study_name, target_vcf, target_tbi, msav, bcf] -> [study_name, chr, target_vcf, target_tbi, ref_bcf, ref_bcf_index]
        ch_phased_input = ch_target_vcf_flat.autosomes.map { chr, study, target, target_tbi, msav, bcf ->
            tuple(
                study,
                chr,
                target,
                target_tbi,
                bcf.find { it.name.endsWith('.bcf') },
                bcf.find { it.name.endsWith('.csi') }
            )
        }
        ch_phased_vcfs = PHASE_AUTOSOMES(ch_phased_input)
    }

    // 3. Imputation using Minimac4
    if (params.run_imputation) {
        // Re-key phased autosomes on chromosome so the join happens against the
        // reference chromosome key, not the msav path.
        ch_impute_input = ch_phased_vcfs
            .map { study, chr, phased, phased_index -> tuple(chr, study, phased, phased_index) }
            .combine(ch_reference.map { chr, msav, bcf ->
                tuple(chr, msav instanceof List ? msav[0] : msav)
            }, by: 0)
            .map { chr, study, phased, phased_index, msav ->
                tuple(study, chr, phased, phased_index, msav)
            }
        
        ch_imputed_vcfs = IMPUTE_AUTOSOMES(ch_impute_input)
        
        // ChrX split, phase, and impute
        ch_impute_chrX_input = ch_target_vcf_flat.chrX.map { chr, study, target, target_tbi, msavs, bcfs ->
            tuple(study, chr, target, target_tbi, msavs, bcfs)
        }
        ch_imputed_chrX = IMPUTE_CHRX(ch_impute_chrX_input)
    }
}
