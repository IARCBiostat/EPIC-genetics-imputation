#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PREP_DBSNP_CHROM } from './modules/prepare_dbsnp.nf'
include { PREP_CHROM } from './modules/prepare_chrom.nf'
include { SAMPLE_REVIEW } from './modules/sample_review.nf'
include { FINALIZE_STUDY } from './modules/finalize_study.nf'

def chromSortKey(String chrom) {
    chrom == 'X' ? 23 : chrom.toInteger()
}

workflow {

    def included_studies = params.study == 'all' ? [] : params.study.split(',').collect { it.trim() }
    def all_chroms = (1..22).collect { it.toString() } + ['X']
    def expected_chrom_count = all_chroms.size()

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

    ch_prep_input = ch_stage2
        .map { study, chr, vcf, tbi -> tuple(chr, study, vcf, tbi) }
        // Fan out each per-chromosome dbSNP reference to every study VCF on that chromosome.
        .combine(ch_dbsnp_chr, by: 0)
        .map { chr, study, vcf, tbi, dbsnp_vcf, dbsnp_tbi ->
            tuple(study, chr, vcf, tbi, dbsnp_vcf, dbsnp_tbi)
        }

    ch_chrom_review = PREP_CHROM(ch_prep_input)

    ch_sample_review_input = ch_chrom_review
        .map { study, chr, pgen, pvar, psam, stats, id_map ->
            tuple(study, [chr: chr, pgen: pgen, pvar: pvar, psam: psam, stats: stats, id_map: id_map])
        }
        .groupTuple(size: expected_chrom_count)
        .map { study, entries ->
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

    ch_study_review = SAMPLE_REVIEW(ch_sample_review_input)

    FINALIZE_STUDY(ch_study_review.finalize_input)
}
