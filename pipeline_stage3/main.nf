#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PREP_DBSNP_CHROM } from './modules/prepare_dbsnp.nf'
include { FILTER_CHROM; ANNOTATE_CHROM; IMPORT_CHROM; HWE_CHROM } from './modules/prepare_chrom.nf'
include { MERGE_STUDY; SEX_CHECK; PRUNE_AUTOSOMES; KING_QC; HET_PCA_QC; SAMPLE_REVIEW_SUMMARY } from './modules/sample_review.nf'
include { FINALIZE_STUDY } from './modules/finalize_study.nf'

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

    // Use per-study group sizes so downstream sample QC starts as soon as a study's
    // own chromosomes are complete, instead of waiting for all studies' HWE_CHROM tasks.
    ch_stage2_chrom_counts = ch_stage2
        .map { study, chr, vcf, tbi -> tuple(study, chr) }
        .groupTuple()
        .map { study, chroms ->
            tuple(study, chroms.collect { it.toString() }.toSet().size())
        }

    // FILTER_CHROM: apply R2/MAF filter and generate sample ID map — 1 CPU, I/O bound.
    ch_filtered = FILTER_CHROM(
        ch_stage2.map { study, chr, vcf, tbi -> tuple(study, chr, vcf, tbi) }
    )

    // ANNOTATE_CHROM: rsID annotation + ID normalisation — 1 CPU, Python bottleneck.
    // Fan out each per-chromosome dbSNP reference to every filtered BCF on that chromosome.
    // id_map is threaded through as a pass-through so IMPORT_CHROM receives it without a join.
    ch_annotated = ANNOTATE_CHROM(
        ch_filtered
            .map { study, chr, filtered_bcf, id_map -> tuple(chr, study, filtered_bcf, id_map) }
            .combine(ch_dbsnp_chr, by: 0)
            .map { chr, study, filtered_bcf, id_map, dbsnp_vcf, dbsnp_tbi ->
                tuple(study, chr, filtered_bcf, dbsnp_vcf, dbsnp_tbi, id_map)
            }
    )
    // ch_annotated emits: (study, chr, annotated_bcf, variant_id_map, id_map)

    // IMPORT_CHROM: PLINK2 BCF → PGEN — 4 CPUs, compute bound.
    // variant_id_map is threaded through as a pass-through so HWE_CHROM receives it without a join.
    ch_imported = IMPORT_CHROM(ch_annotated)
    // ch_imported emits: (study, chr, pgen, pvar, psam, variant_id_map)

    ch_chrom_review = HWE_CHROM(ch_imported)

    ch_sample_review_input = ch_chrom_review
        .map { study, chr, pgen, pvar, psam, stats, id_map ->
            tuple(study, [chr: chr, pgen: pgen, pvar: pvar, psam: psam, stats: stats, id_map: id_map])
        }
        .combine(ch_stage2_chrom_counts, by: 0)
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
                ordered.collect { it.chr },
                ordered.collect { it.pgen },
                ordered.collect { it.pvar },
                ordered.collect { it.psam },
                ordered.collect { it.stats },
                ordered.collect { it.id_map }
            )
        }
        .join(ch_stage1_fam)
        .map { study, chroms, pgens, pvars, psams, stats, id_maps, fam ->
            tuple(study, chroms, pgens, pvars, psams, stats, id_maps, fam)
        }

    ch_merged_study = MERGE_STUDY(ch_sample_review_input).merged

    ch_sexcheck = SEX_CHECK(
        ch_merged_study.map { study, chroms, pgen, pvar, psam, sex_update ->
            tuple(study, chroms, pgen, pvar, psam, sex_update)
        }
    ).sexcheck

    ch_pruned_bed = PRUNE_AUTOSOMES(
        ch_merged_study.map { study, chroms, pgen, pvar, psam, sex_update ->
            tuple(study, pgen, pvar, psam)
        }
    ).pruned_bed

    ch_king = KING_QC(ch_pruned_bed)
    ch_het_pca = HET_PCA_QC(ch_pruned_bed)

    ch_review_summary_input = ch_sexcheck
        .join(ch_king.related_ids)
        .join(ch_het_pca.qc_files)
        .join(
            ch_merged_study.map { study, chroms, pgen, pvar, psam, sex_update ->
                tuple(study, psam)
            }
        )
        .map { study, sexcheck, related_ids, het, eigenvec, eigenval, psam ->
            tuple(study, sexcheck, related_ids, het, eigenvec, eigenval, psam)
        }

    ch_review_summary = SAMPLE_REVIEW_SUMMARY(ch_review_summary_input).summary

    ch_finalize_input = ch_merged_study
        .map { study, chroms, pgen, pvar, psam, sex_update ->
            tuple(study, pgen, pvar, psam, sex_update)
        }
        .join(ch_review_summary)
        .map { study, pgen, pvar, psam, sex_update, remove_list, sample_review_tsv ->
            tuple(study, pgen, pvar, psam, remove_list, sex_update)
        }

    FINALIZE_STUDY(ch_finalize_input)
}
