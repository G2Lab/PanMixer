"""
Pipeline runner: submits all experiment steps sequentially via SLURM dependencies.

Usage:
    python3 run_main.py [--exp EXP] [--fixed_param PARAM] [--strategy STRATEGY]
                        [--capacity_file FILE] [--subjects_file FILE]

Steps:
    1. experiment_starter  (runs locally)
    2. optimize            (sbatch)
    3. stacker             (sbatch, depends on 2)
    4. convert_2_vcf       (sbatch, depends on 3)
    5. combine_vcfs        (sbatch, depends on 4)
    6. gap_score_all       (sbatch, depends on 5)
    7. af_loss             (sbatch, depends on 5)
    8. ld_loss             (sbatch, depends on 5)
    9. vg_prep             (sbatch, depends on 5)
   10. quick_align         (sbatch, depends on 9)
"""

import argparse

from tools.common.experiment_starter import experiment_starter
from tools.common.utils import latest_experiment_number
from tools.panmixer.optimizer import optimizer
from tools.panmixer.stacker import stacker
from tools.common.convert_2_vcf import convert_2_vcf
from tools.common.combine_vcfs import combine_vcfs
from tools.downstream.privacy.diploid_gap_score import diploid_gap_score_computer
from tools.downstream.privacy.MIA_privacy import MIA_privacy_computer
from tools.downstream.utility_in.af_loss import af_loss_computer
from tools.downstream.utility_in.ld_loss import ld_loss
from tools.downstream.utility_out.vg_prep import vg_prep
from tools.downstream.utility_out.quick_align import quick_align
from tools.common.slurm_helper import set_dependency, clear_dependency

from constants import DEFAULT_CAPACITY_FILE, DEFAULT_SUBJECTS_FILE


STEP_ORDER = [
    "optimize", "stacker", "convert_2_vcf", "combine_vcfs",
    "gap_score", "af_loss", "ld_loss", "mia_privacy", "vg_prep", "quick_align",
]


def main():
    parser = argparse.ArgumentParser(description="PanMixer pipeline runner")
    parser.add_argument("--exp", type=int, default=-1, help="Experiment number (-1 = create new)")
    parser.add_argument("--fixed_param", type=str, default="utility", help="Optimizer fixed parameter")
    parser.add_argument("--strategy", type=str, default="to_best", help="Stacker strategy")
    parser.add_argument("--capacity_file", type=str, default=DEFAULT_CAPACITY_FILE)
    parser.add_argument("--subjects_file", type=str, default=DEFAULT_SUBJECTS_FILE)
    parser.add_argument("--seed", type=int, default=None, help="Random seed for sampling steps")
    parser.add_argument("--skip_starter", action="store_true", help="Skip experiment_starter (use existing experiment)")
    parser.add_argument("--stop_after", type=str, default=None, choices=STEP_ORDER,
                        help="Stop submitting steps after this one (inclusive)")
    args = parser.parse_args()

    def should_stop(step_name):
        return args.stop_after is not None and STEP_ORDER.index(step_name) >= STEP_ORDER.index(args.stop_after)

    # --- Step 1: experiment_starter (runs locally, no sbatch) ---
    if not args.skip_starter:
        print("=== Step 1: experiment_starter ===")
        experiment_starter(args.capacity_file, args.subjects_file, False, False, False)

    exp = args.exp
    if exp == -1:
        exp = latest_experiment_number()
    print(f"Using experiment: {exp}")

    submitted = {}

    # --- Step 2: optimize ---
    print("=== Step 2: optimize ===")
    clear_dependency()
    job_optimize = optimizer(False, args.fixed_param, exp, args.seed)
    submitted["optimize"] = job_optimize
    if should_stop("optimize"): return _summary(exp, submitted)

    # --- Step 3: stacker ---
    print("=== Step 3: stacker ===")
    set_dependency(job_optimize)
    job_stacker = stacker(args.strategy, exp)
    submitted["stacker"] = job_stacker
    if should_stop("stacker"): return _summary(exp, submitted)

    # --- Step 4: convert_2_vcf ---
    print("=== Step 4: convert_2_vcf ===")
    set_dependency(job_stacker)
    job_convert = convert_2_vcf(exp)
    submitted["convert_2_vcf"] = job_convert
    if should_stop("convert_2_vcf"): return _summary(exp, submitted)

    # --- Step 5: combine_vcfs ---
    print("=== Step 5: combine_vcfs ===")
    set_dependency(job_convert)
    job_combine = combine_vcfs(exp)
    submitted["combine_vcfs"] = job_combine
    if should_stop("combine_vcfs"): return _summary(exp, submitted)

    # --- Steps 6-8: downstream analyses (parallel, all depend on combine_vcfs) ---
    print("=== Step 6: gap_score_all ===")
    set_dependency(job_combine)
    job_gap = diploid_gap_score_computer(exp)
    submitted["gap_score"] = job_gap
    if should_stop("gap_score"): return _summary(exp, submitted)

    print("=== Step 7: af_loss ===")
    set_dependency(job_combine)
    job_af = af_loss_computer(exp, False)
    submitted["af_loss"] = job_af
    if should_stop("af_loss"): return _summary(exp, submitted)

    print("=== Step 8: ld_loss ===")
    set_dependency(job_combine)
    job_ld = ld_loss(exp, False)
    submitted["ld_loss"] = job_ld
    if should_stop("ld_loss"): return _summary(exp, submitted)

    # --- Step 9: MIA_privacy (depends on stacker, not combine_vcfs) ---
    print("=== Step 9: MIA_privacy ===")
    set_dependency(job_stacker)
    job_mia = MIA_privacy_computer(exp)
    submitted["mia_privacy"] = job_mia
    if should_stop("mia_privacy"): return _summary(exp, submitted)

    # --- Step 10: vg_prep (depends on combine_vcfs) ---
    print("=== Step 10: vg_prep ===")
    set_dependency(job_combine)
    job_vg = vg_prep(exp)
    submitted["vg_prep"] = job_vg
    if should_stop("vg_prep"): return _summary(exp, submitted)

    # --- Step 11: quick_align (depends on vg_prep) ---
    print("=== Step 11: quick_align ===")
    set_dependency(job_vg)
    job_align = quick_align(exp)
    submitted["quick_align"] = job_align

    return _summary(exp, submitted)


def _summary(exp, submitted):
    clear_dependency()
    print("\n=== Submitted jobs ===")
    print(f"Experiment: {exp}")
    for name, jobid in submitted.items():
        print(f"  {name:14s} {jobid}")


if __name__ == "__main__":
    main()
