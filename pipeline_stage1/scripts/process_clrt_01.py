#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))

from script_utils import (
    copy_prefix,
    ensure_dirs,
    liftover_to_hg38,
    manifest_arg,
    move_prefix,
    plink_cmd,
    python2_cmd,
    q,
    remove_prefix,
    run,
    write_empty_file,
    write_summary,
)


PREFIX = 'Clrt_01'
STUDY_ID = 'Clrt_01'
RAW_REL = 'Colonrectum/Clrt_01_Gecco/Data_Received/Data_Extracted/clrt_gecco_geno'
MANIFEST_REL = 'Colonrectum/Clrt_01_Gecco/Chip_files/HumanOmniExpressExome-8-v1-2-B.csv'
MANIFEST_FLAGS = '-c 9 -p 10 -n 1 -a 3 -sp 20 -sf 0 -st 2 -b 8 -nlh 8 -nlt 24'
BUILD = '37'
ID_LINK_REL = 'Colonrectum/Clrt_01_Gecco/Data_Received/EPIC_samplefile.tsv'
ID_FLAGS = '-s 1 -e 2 -nlh 1 -q'
QC_REFERENCE_PREFIX = 'Clrt_01_linked'
PRE = {'name_linkage': True, 'second_qc_mode': 'a'}
PART1 = {'neg': 'exclude_PM', 'chrpos': 'exclude', 'alleles': 'exclude'}
COMPLETION = {'chr': 'move', 'pos': 'move', 'strand': 'move', 'alleles': 'move'}
PART2_BUILD35_REL = None


def render_flags(flags: str, data_root: Path) -> str:
    return flags.replace('{data_root}', str(data_root))


def part1_exclude(plink: str, trace: Path, source_name: str, qcinit: Path, list_name: str, out_name: str, study_dir: Path) -> None:
    run(
        plink_cmd(
            plink,
            f"--bfile {q(trace / source_name)} --exclude {q(qcinit / list_name)} --make-bed --out {q(trace / out_name)}",
        ),
        study_dir,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=f"Stage 1 processing for {STUDY_ID}")
    parser.add_argument('--data-root', required=True, help='Root directory containing the raw EPIC genetics datasets.')
    parser.add_argument('--work-root', default=str(Path(__file__).resolve().parents[1] / 'work'), help='Output root for stage-1 processing results.')
    parser.add_argument('--plink', default='plink', help='PLINK executable to use.')
    parser.add_argument('--python2', default='python2.7', help='Python 2 interpreter for legacy preprocessing utilities.')
    args = parser.parse_args()

    data_root = Path(args.data_root).resolve()
    work_root = Path(args.work_root).resolve()
    pipeline_root = Path(__file__).resolve().parents[1]
    tools = pipeline_root / 'bin'
    study_dir = work_root / STUDY_ID
    trace = study_dir / 'Trace'
    qcinit = trace / 'QC_init'
    qc_completion = trace / 'QC_completion'
    qc_exclusion = trace / 'QC_exclusion'
    lift_dir = trace / 'LiftOver'
    stage1_dir = pipeline_root.parent / 'analysis' / STUDY_ID / 'stage1'
    summary = stage1_dir / 'summary.txt'

    ensure_dirs([study_dir, trace, qcinit, qc_completion, qc_exclusion, lift_dir, stage1_dir])

    manifest_read = manifest_arg(data_root, MANIFEST_REL, render_flags(MANIFEST_FLAGS, data_root))
    base_prefix = trace / PREFIX
    copy_prefix(data_root / RAW_REL, base_prefix, study_dir)

    if PRE.get('sort_input'):
        sort_prefix = trace / f"{PREFIX}_sort"
        run(plink_cmd(args.plink, f"--bfile {q(base_prefix)} --make-bed --out {q(sort_prefix)}"), study_dir)
        remove_prefix(base_prefix, study_dir)
        move_prefix(sort_prefix, base_prefix, study_dir)

    write_summary(summary, "\n*************************** QC of input data *******************************\n", 'w')
    run(
        python2_cmd(args.python2, tools / 'preProcessing.py', f"-sum_n 1 -f {q(base_prefix)} -m {manifest_read} -sum_f {q(summary)}"),
        study_dir,
    )
    run(f"mv {PREFIX}_* {q(qcinit)}", study_dir)
    run('rm -f Ref* Bim', study_dir)

    current_id_prefix_name = PREFIX

    if PRE.get('name_linkage'):
        linked_name = f"{PREFIX}_linked"
        run(
            python2_cmd(args.python2, tools / 'name_linkage_Bim_Manifest.py', f"-f {q(base_prefix)} -m {q(qcinit / (PREFIX + '_manifest.csv'))}"),
            study_dir,
        )
        run(f"mv {PREFIX}_linkManifestNames.txt {q(trace / (PREFIX + '_linkManifestNames.txt'))}", study_dir)
        run(
            plink_cmd(args.plink, f"--bfile {q(base_prefix)} --update-name {q(trace / (PREFIX + '_linkManifestNames.txt'))} --make-bed --out {q(trace / linked_name)}"),
            study_dir,
        )
        run(f"rm -f {q(str(qcinit / PREFIX))}_*", study_dir)
        write_summary(summary, "\n*************************** QC of input data after linkage *******************************\n", PRE.get('second_qc_mode', 'a'))
        run(
            python2_cmd(args.python2, tools / 'preProcessing.py', f"-sum_n 1 -f {q(trace / linked_name)} -m {manifest_read} -sum_f {q(summary)}"),
            study_dir,
        )
        run(f"mv {linked_name}_* {q(qcinit)}", study_dir)
        run('rm -f Ref* Bim', study_dir)
        remove_prefix(base_prefix, study_dir)
        current_id_prefix_name = linked_name

    if PRE.get('reset_positions'):
        run(python2_cmd(args.python2, tools / 'reset_SNP_position.py', f"-f {q(base_prefix)}"), study_dir)
        remove_prefix(base_prefix, study_dir)
        move_prefix(trace / f"{PREFIX}_reset_all_pos", base_prefix, study_dir)
        run(f"rm -f {q(str(qcinit / PREFIX))}_*", study_dir)
        write_summary(summary, "\n*************************** QC of input data after reset *******************************\n", PRE.get('second_qc_mode', 'a'))
        run(
            python2_cmd(args.python2, tools / 'preProcessing.py', f"-sum_n 1 -f {q(base_prefix)} -m {manifest_read} -sum_f {q(summary)}"),
            study_dir,
        )
        run(f"mv {PREFIX}_* {q(qcinit)}", study_dir)
        run('rm -f Ref* Bim', study_dir)

    if PRE.get('exclude_xymt'):
        excxy_name = f"{PREFIX}_EXCXY"
        run(plink_cmd(args.plink, f"--bfile {q(base_prefix)} --not-chr x y xy mt --make-bed --out {q(trace / excxy_name)}"), study_dir)
        run(f"rm -f {q(str(qcinit / PREFIX))}_*", study_dir)
        write_summary(summary, "\n*************************** QC of input data after chromosome exclusion *******************************\n", PRE.get('second_qc_mode', 'a'))
        run(
            python2_cmd(args.python2, tools / 'preProcessing.py', f"-sum_n 1 -f {q(trace / excxy_name)} -m {manifest_read} -sum_f {q(summary)}"),
            study_dir,
        )
        run(f"mv {excxy_name}_* {q(qcinit)}", study_dir)
        run('rm -f Ref* Bim', study_dir)
        remove_prefix(base_prefix, study_dir)
        current_id_prefix_name = excxy_name

    if PRE.get('ab_translate'):
        ab_args = f"-f {q(base_prefix)} -m {q(qcinit / (PREFIX + '_manifest.csv'))}"
        if PRE.get('ab_txt_rel'):
            ab_args += f" -t {q(data_root / PRE['ab_txt_rel'])}"
        run(python2_cmd(args.python2, tools / 'allele_nomenclatureAB_conversion.py', ab_args), study_dir)
        run(f"mv {PREFIX}_alleles_tradAB.txt {q(trace / (PREFIX + '_alleles_tradAB.txt'))}", study_dir)
        trad_ab_name = f"{PREFIX}_tradAB"
        run(
            plink_cmd(args.plink, f"--bfile {q(base_prefix)} --update-alleles {q(trace / (PREFIX + '_alleles_tradAB.txt'))} --make-bed --out {q(trace / trad_ab_name)}"),
            study_dir,
        )
        remove_prefix(base_prefix, study_dir)
        current_id_prefix_name = trad_ab_name

    id_input_prefix = trace / current_id_prefix_name
    id_args = f"-f {q(id_input_prefix.with_suffix('.fam'))} -i {q(data_root / 'Reference/Epic/Subj_Id_2015.txt')}"
    if ID_LINK_REL:
        id_args += f" -l {q(data_root / ID_LINK_REL)}"
    if ID_FLAGS:
        id_args += f" {ID_FLAGS}"
    id_args += f" -sum_f {q(summary)}"

    run(python2_cmd(args.python2, tools / 'link_ID_standard.py', id_args), study_dir)
    run(f"mv {current_id_prefix_name}_goodID.txt {q(trace / (PREFIX + '_goodID.txt'))}", study_dir)
    run(f"mv {current_id_prefix_name}_removeID.txt {q(trace / (PREFIX + '_removeID.txt'))}", study_dir)

    good_ids = trace / f"{PREFIX}_goodID.txt"
    remove_ids = trace / f"{PREFIX}_removeID.txt"
    epic_sex = trace / f"{PREFIX}_EPIC_sex"
    epic_rm = trace / f"{PREFIX}_EPIC_sex_rmID"
    epic_id = trace / f"{PREFIX}_EPIC_sex_ID"
    epic_pheno = trace / f"{PREFIX}_EPIC_sex_ID_pheno"
    epic_all = trace / f"{PREFIX}_EPIC_all"

    run(plink_cmd(args.plink, f"--bfile {q(id_input_prefix)} --update-sex {q(good_ids)} 3 --make-bed --out {q(epic_sex)}"), study_dir)
    run(plink_cmd(args.plink, f"--bfile {q(epic_sex)} --remove {q(remove_ids)} --make-bed --out {q(epic_rm)}"), study_dir)
    run(plink_cmd(args.plink, f"--bfile {q(epic_rm)} --update-ids {q(good_ids)} --make-bed --out {q(epic_id)}"), study_dir)
    empty_file = trace / 'empty_file.txt'
    write_empty_file(empty_file)
    run(plink_cmd(args.plink, f"--bfile {q(epic_id)} --make-pheno {q(empty_file)} 1 --make-bed --out {q(epic_pheno)}"), study_dir)
    empty_file.unlink(missing_ok=True)
    run(plink_cmd(args.plink, f"--bfile {q(epic_pheno)} --indiv-sort n --make-bed --out {q(epic_all)}"), study_dir)
    run(f"rm -f {q(str(trace / (PREFIX + '_EPIC_sex')))}*", study_dir)
    remove_prefix(id_input_prefix, study_dir)

    part1_current = f"{PREFIX}_EXC1_SNPfound"
    part1_exclude(args.plink, trace, f"{PREFIX}_EPIC_all", qcinit, f"{QC_REFERENCE_PREFIX}_not_found.txt", part1_current, study_dir)

    neg_action = PART1.get('neg')
    if neg_action == 'move':
        next_name = f"{PREFIX}_EXC1_SNPfound_strand"
        move_prefix(trace / part1_current, trace / next_name, study_dir)
        part1_current = next_name
    elif neg_action and neg_action.startswith('exclude_'):
        next_name = f"{PREFIX}_EXC1_SNPfound_strand"
        strand_tag = neg_action.replace('exclude_', '')
        part1_exclude(args.plink, trace, part1_current, qcinit, f"{QC_REFERENCE_PREFIX}_negative_strand_{strand_tag}.txt", next_name, study_dir)
        part1_current = next_name

    chrpos_action = PART1.get('chrpos')
    if chrpos_action == 'move':
        next_name = f"{PREFIX}_EXC1_SNPfound_strand_samePos"
        move_prefix(trace / part1_current, trace / next_name, study_dir)
        part1_current = next_name
    elif chrpos_action == 'exclude':
        next_name = f"{PREFIX}_EXC1_SNPfound_strand_samePos"
        part1_exclude(args.plink, trace, part1_current, qcinit, f"{QC_REFERENCE_PREFIX}_diff_chrpos.txt", next_name, study_dir)
        part1_current = next_name

    allele_action = PART1.get('alleles')
    if allele_action == 'move_all':
        next_name = f"{PREFIX}_EXC1_all"
        move_prefix(trace / part1_current, trace / next_name, study_dir)
        part1_current = next_name
    elif allele_action == 'move_same':
        next_name = f"{PREFIX}_EXC1_SNPfound_strand_samePos_sameAllele"
        move_prefix(trace / part1_current, trace / next_name, study_dir)
        part1_current = next_name
    elif allele_action == 'exclude':
        next_name = f"{PREFIX}_EXC1_all"
        part1_exclude(args.plink, trace, part1_current, qcinit, f"{QC_REFERENCE_PREFIX}_diff_alleles.txt", next_name, study_dir)
        part1_current = next_name
    elif allele_action == 'exclude_same':
        next_name = f"{PREFIX}_EXC1_SNPfound_strand_samePos_sameAllele"
        part1_exclude(args.plink, trace, part1_current, qcinit, f"{QC_REFERENCE_PREFIX}_diff_alleles.txt", next_name, study_dir)
        part1_current = next_name

    if PART1.get('unknown'):
        next_name = f"{PREFIX}_EXC1_all"
        part1_exclude(args.plink, trace, part1_current, qcinit, f"{QC_REFERENCE_PREFIX}_unknown_strand_{PART1['unknown']}.txt", next_name, study_dir)
        part1_current = next_name

    run(f"rm -f {q(str(trace / (PREFIX + '_EXC1_SNPfound')))}*", study_dir)
    remove_prefix(epic_all, study_dir)

    completion_manifest = qcinit / f"{QC_REFERENCE_PREFIX}_manifest.csv"
    completion_current = part1_current
    completed_chr = f"{PREFIX}_completed_chr"
    completed_pos = f"{PREFIX}_completed_chr_pos"
    completed_strand = f"{PREFIX}_completed_chr_pos_strand"
    completed_all = f"{PREFIX}_completed_all"
    completion_txt_arg = ''
    if COMPLETION.get('txt_rel'):
        completion_txt_arg = f" -t {q(data_root / COMPLETION['txt_rel'])}"

    if COMPLETION.get('chr') == 'search':
        run(python2_cmd(args.python2, tools / 'search_chr_standard.py', f"-b {q((trace / completion_current).with_suffix('.bim'))} -m {q(completion_manifest)}{completion_txt_arg} -bd {BUILD}"), study_dir)
        run(f"mv {completion_current}_goodChr_{BUILD}.txt {q(trace / (PREFIX + '_completed_goodChr_' + BUILD + '.txt'))}", study_dir)
        run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --update-chr {q(trace / (PREFIX + '_completed_goodChr_' + BUILD + '.txt'))} --make-bed --out {q(trace / completed_chr)}"), study_dir)
        remove_prefix(trace / completion_current, study_dir)
    else:
        move_prefix(trace / completion_current, trace / completed_chr, study_dir)
    completion_current = completed_chr

    if COMPLETION.get('pos') == 'search':
        run(python2_cmd(args.python2, tools / 'search_pos_standard.py', f"-b {q((trace / completion_current).with_suffix('.bim'))} -m {q(completion_manifest)}{completion_txt_arg} -bd {BUILD}"), study_dir)
        run(f"mv {completion_current}_goodPos_{BUILD}.txt {q(trace / (completion_current + '_goodPos_' + BUILD + '.txt'))}", study_dir)
        run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --update-map {q(trace / (completion_current + '_goodPos_' + BUILD + '.txt'))} --make-bed --out {q(trace / completed_pos)}"), study_dir)
        remove_prefix(trace / completion_current, study_dir)
    else:
        move_prefix(trace / completion_current, trace / completed_pos, study_dir)
    completion_current = completed_pos

    strand_mode = COMPLETION.get('strand')
    if strand_mode == 'flip_tb_pm':
        temp_name = f"{PREFIX}_completed_chr_pos_strand1"
        run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --flip {q(qcinit / (QC_REFERENCE_PREFIX + '_negative_strand_manifest_TB.txt'))} --make-bed --out {q(trace / temp_name)}"), study_dir)
        run(plink_cmd(args.plink, f"--bfile {q(trace / temp_name)} --flip {q(qcinit / (QC_REFERENCE_PREFIX + '_negative_strand_manifest_PM.txt'))} --make-bed --out {q(trace / completed_strand)}"), study_dir)
        remove_prefix(trace / temp_name, study_dir)
        remove_prefix(trace / completion_current, study_dir)
    elif strand_mode == 'flip_fr_pm':
        temp_name = f"{PREFIX}_completed_chr_pos_strand1"
        run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --flip {q(qcinit / (QC_REFERENCE_PREFIX + '_negative_strand_manifest_FR.txt'))} --make-bed --out {q(trace / temp_name)}"), study_dir)
        run(plink_cmd(args.plink, f"--bfile {q(trace / temp_name)} --flip {q(qcinit / (QC_REFERENCE_PREFIX + '_negative_strand_manifest_PM.txt'))} --make-bed --out {q(trace / completed_strand)}"), study_dir)
        remove_prefix(trace / temp_name, study_dir)
        remove_prefix(trace / completion_current, study_dir)
    elif strand_mode == 'flip_pm':
        run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --flip {q(qcinit / (QC_REFERENCE_PREFIX + '_negative_strand_manifest_PM.txt'))} --make-bed --out {q(trace / completed_strand)}"), study_dir)
        remove_prefix(trace / completion_current, study_dir)
    else:
        move_prefix(trace / completion_current, trace / completed_strand, study_dir)
    completion_current = completed_strand

    if COMPLETION.get('alleles') == 'search':
        run(python2_cmd(args.python2, tools / 'search_alleles_standard.py', f"-b {q((trace / completion_current).with_suffix('.bim'))} -m {q(completion_manifest)}{completion_txt_arg} -bd {BUILD}"), study_dir)
        run(f"mv {completion_current}_goodAlleles_{BUILD}.txt {q(trace / (completion_current + '_goodAlleles_' + BUILD + '.txt'))}", study_dir)
        if COMPLETION.get('subversion_rel'):
            interim = f"{PREFIX}_completed_chr_pos_strand_allel"
            run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --update-alleles {q(trace / (completion_current + '_goodAlleles_' + BUILD + '.txt'))} --make-bed --out {q(trace / interim)}"), study_dir)
            remove_prefix(trace / completion_current, study_dir)
            run(python2_cmd(args.python2, tools / 'subVersion_position_correction.py', f"-b {q((trace / interim).with_suffix('.bim'))} -t {q(qcinit / COMPLETION['subversion_rel'])} -d 1"), study_dir)
            sub_base = COMPLETION['subversion_rel'].replace('.txt', '')
            run(f"mv {sub_base}_goodSubversionPos.txt {q(trace / (sub_base + '_goodSubversionPos.txt'))}", study_dir)
            run(plink_cmd(args.plink, f"--bfile {q(trace / interim)} --update-map {q(trace / (sub_base + '_goodSubversionPos.txt'))} --make-bed --out {q(trace / completed_all)}"), study_dir)
            remove_prefix(trace / interim, study_dir)
        else:
            run(plink_cmd(args.plink, f"--bfile {q(trace / completion_current)} --update-alleles {q(trace / (completion_current + '_goodAlleles_' + BUILD + '.txt'))} --make-bed --out {q(trace / completed_all)}"), study_dir)
            remove_prefix(trace / completion_current, study_dir)
    else:
        move_prefix(trace / completion_current, trace / completed_all, study_dir)

    write_summary(summary, "\n*************************** QC after completing data *******************************\n", 'a')
    run(python2_cmd(args.python2, tools / 'preProcessing.py', f"-sum_n 2 -f {q(trace / completed_all)} -m {manifest_read} -sum_f {q(summary)}"), study_dir)
    run(f"mv {completed_all}_* {q(qc_completion)}", study_dir)
    run('rm -f Ref* Bim', study_dir)
    run(f"rm -f {q(str(qc_completion / (completed_all + '_manifest.csv')))}", study_dir)
    run(f"rm -f {q(str(qc_completion / (completed_all + '_negative_strand_manifest_')))}*", study_dir)

    exc2_mito = f"{PREFIX}_EXC2_mito"
    exc2_t4 = f"{PREFIX}_EXC2_mito_T4"
    exc2_t2 = f"{PREFIX}_EXC2_mito_T4_T2"
    exc2_t3 = f"{PREFIX}_EXC2_mito_T4_T2_T3"
    final_exc2 = f"{PREFIX}_EXC2_all"

    part1_exclude(args.plink, trace, completed_all, qc_completion, f"{completed_all}_missplaced_mito.txt", exc2_mito, study_dir)
    part1_exclude(args.plink, trace, exc2_mito, qc_completion, f"{completed_all}_duplicate_T4_error.txt", exc2_t4, study_dir)
    part1_exclude(args.plink, trace, exc2_t4, qc_completion, f"{completed_all}_duplicate_T2_name_pos.txt", exc2_t2, study_dir)
    part1_exclude(args.plink, trace, exc2_t2, qc_completion, f"{completed_all}_duplicate_T3_pos.txt", exc2_t3, study_dir)

    if PART2_BUILD35_REL:
        t1_out = f"{PREFIX}_EXC2_prebuild35"
        part1_exclude(args.plink, trace, exc2_t3, qc_completion, f"{completed_all}_duplicate_T1_same_allele.txt", t1_out, study_dir)
        part1_exclude(args.plink, trace, t1_out, qcinit, PART2_BUILD35_REL, final_exc2, study_dir)
        remove_prefix(trace / t1_out, study_dir)
    else:
        part1_exclude(args.plink, trace, exc2_t3, qc_completion, f"{completed_all}_duplicate_T1_same_allele.txt", final_exc2, study_dir)

    run(f"rm -f {q(str(trace / (PREFIX + '_EXC2_mito')))}*", study_dir)

    run(f"find {q(qc_completion)} -type f -empty -delete", study_dir)

    write_summary(summary, "\n*************************** QC after Excluding SNP Part2 *******************************\n", 'a')
    run(python2_cmd(args.python2, tools / 'preProcessing.py', f"-sum_n 3 -f {q(trace / final_exc2)} -m {manifest_read} -sum_f {q(summary)}"), study_dir)
    run(f"mv {final_exc2}_* {q(qc_exclusion)}", study_dir)
    run('rm -f Ref* Bim', study_dir)
    run(f"rm -f {q(str(qc_exclusion / (final_exc2 + '_manifest.csv')))}", study_dir)
    run(f"rm -f {q(str(qc_exclusion / (final_exc2 + '_negative_strand_manifest_')))}*", study_dir)
    run(f"find {q(qc_exclusion)} -type f -empty -delete", study_dir)
    remove_prefix(trace / completed_all, study_dir)

    liftover_to_hg38(
        args.plink,
        pipeline_root.parent,
        trace / final_exc2,
        stage1_dir / STUDY_ID,
        BUILD,
        lift_dir,
        summary,
        study_dir,
    )


if __name__ == '__main__':
    main()
