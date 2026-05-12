#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * MAIN EPIC IMPUTATION PIPELINE
 * Input datasets come from stage 1 outputs in analysis/<STUDY>/stage1/
 */

include { PREP_REFERENCE } from './modules/reference_prep.nf'
include { PREP_TARGET_VCF } from './modules/target_prep.nf'
include { PHASE_AUTOSOMES } from './modules/phasing.nf'
include { STAGE_PHASING_METRICS } from './modules/phasing_metrics.nf'
include { IMPUTE_AUTOSOMES } from './modules/imputation.nf'
include { IMPUTE_CHRX } from './modules/chrX_imputation.nf'
include { STAGE_IMPUTATION_METRICS } from './modules/summary_metrics.nf'
include { STAGE_EMPIRICAL_VALIDATION_METRICS } from './modules/empirical_validation.nf'
include { REPORTING } from './modules/reporting.nf'

workflow {
    ch_reporting_inputs = Channel.empty()

    // Define requested chromosomes. Use --chromosomes 22 for a single-chromosome smoke run.
    def all_chroms = (1..22).collect { it.toString() } + ['X']
    def selected_chroms = params.chromosomes == 'all'
        ? all_chroms
        : params.chromosomes.toString().tokenize(',').collect { it.trim().replaceFirst(/^chr/, '') }.collect { it == '23' ? 'X' : it }
    selected_chroms.each { chrom ->
        if (!all_chroms.contains(chrom)) {
            error "Unsupported chromosome '${chrom}'. Use 1..22, X, chrX, 23, or all."
        }
    }

    ch_chroms_auto = Channel.fromList(selected_chroms.findAll { it != 'X' })
    ch_chroms_X = Channel.fromList(selected_chroms.findAll { it == 'X' })
    def run_empirical_validation = params.run_empirical_validation.toString().toBoolean()

    // Parse included studies
    def included_studies = params.study == 'all' ? [] : params.study.split(',').collect{ it.trim() }

    def stage1_glob = "${params.stage1_root}/*/stage1/*.{bed,bim,fam}"

    // Read authoritative final stage-1 handoff files.
    ch_stage1_plink = Channel.fromFilePairs(stage1_glob, size: 3)
        .map { prefix, files ->
            def bed = files.find { it.name.endsWith('.bed') }
            def bim = files.find { it.name.endsWith('.bim') }
            def fam = files.find { it.name.endsWith('.fam') }
            def studyName = bed.baseName
            tuple(studyName, prefix, bed, bim, fam)
        }
        .filter { study_name, prefix, bed, bim, fam -> 
            params.study == 'all' || included_studies.contains(study_name) 
        }

    // 1. Reference pre-processing and stage-1 target export
    if (params.run_preprocessing) {
        ch_reference = PREP_REFERENCE(ch_chroms_auto.mix(ch_chroms_X))

        // Convert stage-1 hg38 PLINK inputs to one indexed target VCF per study/chromosome.
        // Chromosome is an explicit input so a chr22 smoke run cannot satisfy an all-chromosome rerun from cache.
        ch_target_input = ch_stage1_plink
            .combine(ch_chroms_auto.mix(ch_chroms_X))
            .map { study_name, prefix, bed, bim, fam, chr ->
                tuple(study_name, chr, prefix, bed, bim, fam)
            }
        ch_target_vcf = PREP_TARGET_VCF(ch_target_input)

        ch_target_vcf_flat = ch_target_vcf
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
        // ch_phasing_metrics = STAGE_PHASING_METRICS(ch_phased_vcfs)
        // ch_reporting_inputs = ch_reporting_inputs.mix(ch_phasing_metrics)
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

        ch_imputed_chrX_for_metrics = ch_imputed_chrX.map { study, vcf, tbi ->
            tuple(study, 'X', vcf, tbi)
        }
        ch_imputation_metrics = STAGE_IMPUTATION_METRICS(
            ch_imputed_vcfs.mix(ch_imputed_chrX_for_metrics)
        )
        ch_reporting_inputs = ch_reporting_inputs.mix(ch_imputation_metrics)

        if (run_empirical_validation) {
            ch_validation_targets = ch_target_vcf.map { chr, study, target, target_tbi ->
                tuple("${study}::${chr}", study, chr, target, target_tbi)
            }
            ch_validation_imputed = ch_imputed_vcfs
                .mix(ch_imputed_chrX_for_metrics)
                .map { study, chr, imputed, imputed_tbi ->
                    tuple("${study}::${chr}", study, chr, imputed, imputed_tbi)
                }
            ch_validation_input = ch_validation_targets
                .combine(ch_validation_imputed, by: 0)
                .map { key, study, chr, target, target_tbi, imputed_study, imputed_chr, imputed, imputed_tbi ->
                    tuple(study, chr, target, target_tbi, imputed, imputed_tbi)
                }
            ch_empirical_validation_metrics = STAGE_EMPIRICAL_VALIDATION_METRICS(ch_validation_input)
            ch_reporting_inputs = ch_reporting_inputs.mix(ch_empirical_validation_metrics)
        }
    }

    if (params.run_reporting) {
        // Each study sends 1 tuple per chromosome from imputation metrics,
        // plus 1 tuple per chromosome from empirical validation if enabled.
        def expected_per_study = selected_chroms.size() * (run_empirical_validation ? 2 : 1)
        
        ch_reporting_inputs
            .groupTuple(by: 0, size: expected_per_study, remainder: true)
            .map { study, metrics, summaries -> tuple(study, metrics + summaries) }
            .set { ch_reporting_grouped }

        REPORTING(ch_reporting_grouped)
    }
}
