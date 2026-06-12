#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PREP_DBSNP_CHROM } from './modules/prepare_dbsnp.nf'
include { PREPARE_CHROM; IMPORT_CHROM; HWE_CHROM; AGGREGATE_HWE_EXCLUDE } from './modules/prepare_chrom.nf'
include { PRUNE_CHROM_FOR_QC; MERGE_PRUNED_FOR_QC; MAKE_UPDATE_FILES; KING_QC; HET_PCA_QC; SAMPLE_REVIEW_SUMMARY } from './modules/sample_review.nf'
include { FINALIZE_CHROM } from './modules/finalize_study.nf'

def chromSortKey(String chrom) {
    chrom == 'X' ? 23 : chrom.toInteger()
}

workflow {

    def included_studies = params.study == 'all' ? [] : params.study.split(',').collect { it.trim() }
    def all_chroms = (1..22).collect { it.toString() } + ['X']

    ch_chroms = Channel.fromList(all_chroms)

    ch_dbsnp_chr = PREP_DBSNP_CHROM(
        ch_chroms.map { chr -> tuple(chr, file(params.dbsnp_vcf), file(params.dbsnp_tbi)) }
    )

    ch_stage2 = Channel.fromPath("${params.stage2_root}/*/stage2/*_chr*_GxS.imputed.vcf.gz")
        .map { vcf ->
            def matcher = (vcf.name =~ /^(.*)_chr([0-9X]+)_GxS\.imputed\.vcf\.gz$/)
            assert matcher.matches() : "Could not parse study/chromosome from ${vcf.name}"
            tuple(matcher[0][1], matcher[0][2], vcf, file("${vcf}.tbi"))
        }
        .filter { study, chr, vcf, tbi ->
            params.study == 'all' || included_studies.contains(study)
        }

    ch_stage1_fam = Channel.fromPath("${params.stage1_root}/*/stage1/*.fam")
        .map { fam -> tuple(fam.baseName, fam) }
        .filter { study, fam ->
            params.study == 'all' || included_studies.contains(study)
        }

    // Fork ch_stage1_fam: consumed by HWE_CHROM (operator) and MAKE_UPDATE_FILES (operator).
    ch_stage1_fam
        .multiMap { study, fam ->
            for_hwe:    tuple(study, fam)
            for_update: tuple(study, fam)
        }
        .set { ch_fam }

    // Per-study total chromosome counts — consumed by AGGREGATE_HWE_EXCLUDE groupTuple.
    ch_stage2
        .map { study, chr, vcf, tbi -> tuple(study, chr) }
        .groupTuple()
        .map { study, chroms ->
            tuple(study, chroms.collect { it.toString() }.toSet().size())
        }
        .set { ch_chrom_counts }

    // Per-study autosome counts (chr 1-22) — consumed by MERGE_PRUNED_FOR_QC groupTuple.
    ch_stage2
        .map { study, chr, vcf, tbi -> tuple(study, chr) }
        .filter { study, chr -> chr != 'X' }
        .groupTuple()
        .map { study, chroms ->
            tuple(study, chroms.collect { it.toString() }.toSet().size())
        }
        .set { ch_autosome_counts }

    // PREPARE_CHROM: R2/MAF filter + rsID annotation in one BCF pass — 4 CPUs.
    ch_prepared = PREPARE_CHROM(
        ch_stage2
            .map { study, chr, vcf, tbi -> tuple(chr, study, vcf, tbi) }
            .combine(ch_dbsnp_chr, by: 0)
            .map { chr, study, vcf, tbi, dbsnp_vcf, dbsnp_tbi ->
                tuple(study, chr, vcf, tbi, dbsnp_vcf, dbsnp_tbi)
            }
    )
    // ch_prepared emits: (study, chr, annotated_bcf, variant_id_map, id_map)

    // IMPORT_CHROM: BCF → PGEN — 8 CPUs, compute bound.
    ch_imported = IMPORT_CHROM(ch_prepared)
    // ch_imported emits: (study, chr, pgen, pvar, psam, variant_id_map)

    // HWE_CHROM: optional HWE exclusion list on controls only — 1 CPU.
    // Fan out the per-study FAM to each chromosome.
    ch_hwe_input = ch_imported
        .combine(ch_fam.for_hwe, by: 0)
    // ch_hwe_input: (study, chr, pgen, pvar, psam, variant_id_map, fam)

    ch_hwe_out = HWE_CHROM(ch_hwe_input)
    ch_chrom_review = ch_hwe_out.chrom_data
    // ch_chrom_review emits: (study, chr, pgen, pvar, psam, variant_metrics)

    // Fork ch_chrom_review for three consumers.
    ch_chrom_review
        .multiMap { study, chr, pgen, pvar, psam, metrics ->
            for_prune:    tuple(study, chr, pgen, pvar, psam)
            for_finalize: tuple(study, chr, pgen, pvar, psam, metrics)
            for_update:   tuple(study, chr, psam)
        }
        .set { ch_chrom_review_split }

    // ── Per-chromosome LD pruning (autosomes only) ────────────────────────────────

    ch_pruned_chrs = PRUNE_CHROM_FOR_QC(
        ch_chrom_review_split.for_prune
            .filter { study, chr, pgen, pvar, psam -> chr != 'X' }
    )

    // Collect 22 per-chromosome pruned BEDs per study, then merge into one small BED.
    ch_merge_pruned_input = ch_pruned_chrs.pruned_chr
        .map { study, chr, bed, bim, fam ->
            tuple(study, [chr: chr, bed: bed, bim: bim, fam: fam])
        }
        .combine(ch_autosome_counts, by: 0)
        .map { study, entry, expected_count ->
            tuple(groupKey(study, expected_count as int), entry)
        }
        .groupTuple()
        .map { study_key, entries ->
            def study = study_key.getGroupTarget()
            def ordered = entries.sort { a, b ->
                chromSortKey(a.chr as String) <=> chromSortKey(b.chr as String)
            }
            tuple(
                study,
                ordered.collect { it.bed },
                ordered.collect { it.bim },
                ordered.collect { it.fam }
            )
        }

    ch_pruned_bed = MERGE_PRUNED_FOR_QC(ch_merge_pruned_input).pruned_bed

    ch_king     = KING_QC(ch_pruned_bed)
    ch_het_pca  = HET_PCA_QC(ch_pruned_bed)

    // ── Sex/pheno update files from chr1 psam + stage1 FAM ───────────────────────
    // All per-chromosome psams are sample-identical; chr1 is the canonical source.

    ch_chr1_psam = ch_chrom_review_split.for_update
        .filter { study, chr, psam -> chr == '1' }
        .map    { study, chr, psam -> tuple(study, psam) }

    // Split chr1_psam for two consumers: MAKE_UPDATE_FILES and SAMPLE_REVIEW_SUMMARY.
    ch_chr1_psam
        .multiMap { study, psam ->
            for_update: tuple(study, psam)
            for_review: tuple(study, psam)
        }
        .set { ch_chr1_psam_split }

    ch_update_files = MAKE_UPDATE_FILES(
        ch_chr1_psam_split.for_update.join(ch_fam.for_update)
    )

    // ── Sample review ─────────────────────────────────────────────────────────────

    ch_review_summary_input = ch_king.related_ids
        .join(ch_het_pca.qc_files)
        .join(ch_chr1_psam_split.for_review)
        .map { study, related_ids, het, eigenvec, eigenval, psam ->
            tuple(study, related_ids, het, eigenvec, eigenval, psam)
        }

    ch_review = SAMPLE_REVIEW_SUMMARY(ch_review_summary_input, params.ancestry, params.related)

    // ── Per-chromosome finalization ───────────────────────────────────────────────
    // Fan out samples_to_remove, sex_update, pheno_update to each FINALIZE_CHROM task.

    ch_finalize_input = ch_chrom_review_split.for_finalize
        .map { study, chr, pgen, pvar, psam, metrics ->
            tuple(study, chr, pgen, pvar, psam)
        }
        .combine(
            ch_review.summary.map { study, remove_list, review_tsv -> tuple(study, remove_list) },
            by: 0
        )
        .combine(
            ch_update_files.update_files.map { study, sex_update, pheno_update -> tuple(study, sex_update, pheno_update) },
            by: 0
        )
        .map { study, chr, pgen, pvar, psam, remove_list, sex_update, pheno_update ->
            tuple(study, chr, pgen, pvar, psam, remove_list, sex_update, pheno_update)
        }

    FINALIZE_CHROM(ch_finalize_input)

    // ── Aggregate per-chromosome HWE exclusion lists into one per-study file ─────

    ch_hwe_study = ch_hwe_out.hwe_exclude
        .combine(ch_chrom_counts, by: 0)
        .map { study, file, expected_count ->
            tuple(groupKey(study, expected_count as int), file)
        }
        .groupTuple()
        .map { study_key, files -> tuple(study_key.getGroupTarget(), files) }

    AGGREGATE_HWE_EXCLUDE(ch_hwe_study)
}
