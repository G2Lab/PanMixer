import numpy as np
import argparse

from tools.common.experiment_starter import experiment_starter
from tools.common.utils import latest_experiment_number
from tools.panmixer.optimizer import optimizer
from tools.panmixer.stacker import stacker
from tools.common.gather_results import gather_results
from tools.downstream.privacy.gap_score import gap_score_computer
from tools.common.convert_2_vcf import convert_2_vcf
from tools.beagle.beagle_refinement import beagle_refinement
from tools.common.combine_vcfs import combine_vcfs
from tools.downstream.utility_in.af_loss import af_loss_computer
from tools.downstream.utility_in.ld_loss import ld_loss
from tools.common.verify_genotypes import verify_all
from tools.beagle.accuracy_stats import accuracy_stats
from tools.downstream.utility_out.vg_prep import vg_prep
from tools.downstream.utility_out.quick_align import quick_align
from tools.downstream.utility_out.filtered_read_mapping import filtered_read_mapping
from tools.downstream.utility_out.personalized_read_mapping import personalized_read_mapping
from tools.downstream.privacy.diploid_gap_score import diploid_gap_score_computer
from tools.downstream.privacy.MIA_privacy import MIA_privacy_computer
from tools.common.create_multitarget_vcfs import create_multitarget_vcfs
from tools.common.VCFtoNP_parallel import VCFtoNP_parallel

from constants import (
    STARTING_DATA_PATH,
    DEFAULT_CAPACITY_FILE,
    DEFAULT_SUBJECTS_FILE,
)

def main():
    parser = argparse.ArgumentParser(description="PanMixer4 Toolkit")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing analysis / files")
    parser.add_argument("--exp" , type=int, default=-1, help="Experiment number")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for tools that sample")

    tools = parser.add_subparsers(dest="tool", required=True, help="Tool to run")
    
    experiment_starter_parser = tools.add_parser("experiment_starter", help="Start an experiment")
    experiment_starter_parser.add_argument("--capacity_file", type=str, default=DEFAULT_CAPACITY_FILE, help="List of capacity constrains")
    experiment_starter_parser.add_argument("--subjects_file", type=str, default=DEFAULT_SUBJECTS_FILE, help="List of subjects")
    experiment_starter_parser.add_argument("--baseline_unedited", action="store_true", help="Run baseline unedit experiment", default=False)
    experiment_starter_parser.add_argument("--baseline_empty", action="store_true", help="Run baseline empty experiment", default=False)
    experiment_starter_parser.add_argument("--baseline_unique", action="store_true", help="Run baseline remove unique experiment", default=False)

    linear_optimizer_parser = tools.add_parser("optimize", help="Run linear optimizer to find optimal obfuscation")
    linear_optimizer_parser.add_argument("--fixed_param", type=str, default="utility", help="Parameter to fix")
    linear_optimizer_parser.add_argument("--baseline_unique", action="store_true", help="Optimize for baseline with removed unique variants", default=False)

    gather_results_parser = tools.add_parser("gather_results", help="Gather results")
    gather_results_parser.add_argument("--optimizer", action="store_true", help="Gather just optimizer results", default=False)
    gather_results_parser.add_argument("--reindex", action="store_true", help="Index experiment", default=False)
    gather_results_parser.add_argument("--gap_score", action="store_true", help="Gather just stacker results", default=False)
    gather_results_parser.add_argument("--stacker", action="store_true", help="Gather just stacker results", default=False)
    gather_results_parser.add_argument("--af_loss", action="store_true", help="Gather just AF loss results", default=False)
    gather_results_parser.add_argument("--ld_loss", action="store_true", help="Gather just AF loss results", default=False)
    gather_results_parser.add_argument("--pangenie_stats", action="store_true", help="Gather just pangenie stats results", default=False)
    gather_results_parser.add_argument("--accuracy_stats", action="store_true", help="Gather just accuracy stats results", default=False)
    gather_results_parser.add_argument("--giraffe", action="store_true", help="Gather just giraffe", default=False)
    gather_results_parser.add_argument("--filtered_giraffe", action="store_true", help="Gather just filtered giraffe", default=False)
    gather_results_parser.add_argument("--personalized_giraffe", action="store_true", help="Gather just personalized giraffe", default=False)
    gather_results_parser.add_argument("--MIA_privacy", action="store_true", help="Gather just MIA privacy results", default=False)

    MIA_privacy_parser = tools.add_parser("MIA_privacy", help="Compute MIA privacy scores")

    stacker_parser = tools.add_parser("stacker", help="Stacker")
    stacker_parser.add_argument("--strategy", type=str, required=True, help="Stacker strategy")

    gapscore_parser = tools.add_parser("gap_score", help="Compute gap score")
    diploid_gapscore_parser = tools.add_parser("gap_score_all", help="Compute gap score")

    af_loss_parser = tools.add_parser("af_loss", help="Compute AF loss")
    af_loss_parser.add_argument("--dont_replace", action="store_true", help="Don't replace", default=False)

    convert_2_vcf_parser = tools.add_parser("convert_2_vcf", help="Convert numpy to vcf")

    beagle_refinement_parser = tools.add_parser("beagle", help="Run beagle refinement")
    
    beagle_stats_parser = tools.add_parser("beagle_stats", help="Run beagle stats")

    ld_loss_parser = tools.add_parser("ld_loss", help="Compute LD loss")
    ld_loss_parser.add_argument("--dont_replace", action="store_true", help="Don't replace", default=False)

    get_alignment_accuracy = tools.add_parser("get_alignment_accuracy", help="Get alignment accuracy")

    accuracy_stats_parser = tools.add_parser("accuracy_stats", help="Get accuracy stats")
    
    combined_vcf_parser = tools.add_parser("combine_vcfs", help="Combine chromosome vcf into one for aligning")

    verify_all_parser = tools.add_parser("verify", help="Verify all samples")

    vg_prep_parser = tools.add_parser("vg_prep", help="Preps giraffe alignment with vg")

    quick_align_parser = tools.add_parser("quick_align", help="Quick alignment with chr21 reads")

    filtered_read_mapping_parser = tools.add_parser("filtered_read_mapping", help="Read mapping to AF>=10%% filtered graph")

    personalized_read_mapping_parser = tools.add_parser("personalized_read_mapping", help="Read mapping to personalized (unfiltered) graph")

    create_multitarget_vcfs_parser = tools.add_parser("create_multitarget_vcfs", help="Create multitarget vcfs")
    create_multitarget_vcfs_parser.add_argument("--target_exp", type=int, required=True, help="Target experiment")

    VCFtoNP_parallel_parser = tools.add_parser("VCFtoNP_parallel", help="Convert vcf to np")

    args = parser.parse_args()

    exp = args.exp
    if exp == -1:
        exp = latest_experiment_number()
    if args.tool == "experiment_starter":
        experiment_starter(args.capacity_file, args.subjects_file, args.baseline_unedited, args.baseline_empty, args.baseline_unique)
    elif args.tool == "optimize":
        optimizer(args.baseline_unique, args.fixed_param, exp, args.seed)
    elif args.tool == "gather_results":
        gather_results(exp, args.overwrite, args.optimizer, args.reindex, args.gap_score, args.stacker, args.af_loss, args.pangenie_stats, args.ld_loss, args.accuracy_stats, args.giraffe, args.filtered_giraffe, args.personalized_giraffe, args.MIA_privacy)
    elif args.tool == "stacker":
        stacker(args.strategy, exp)
    elif args.tool == "gap_score":
        gap_score_computer(exp)
    elif args.tool == "gap_score_all":
        diploid_gap_score_computer(exp)
    elif args.tool == "MIA_privacy":
        MIA_privacy_computer(exp)
    elif args.tool == "convert_2_vcf":
        convert_2_vcf(exp)
    elif args.tool == "beagle":
        beagle_refinement(exp)
    elif args.tool == "combine_vcfs":
        combine_vcfs(exp)
    elif args.tool == "af_loss":
        af_loss_computer(exp, args.dont_replace)
    elif args.tool == "ld_loss":
        ld_loss(exp, args.dont_replace)
    elif args.tool == "verify":
        verify_all(exp)
    elif args.tool == "accuracy_stats":
        accuracy_stats(exp)
    elif args.tool == "vg_prep":
        vg_prep(exp)
    elif args.tool == "quick_align":
        quick_align(exp)
    elif args.tool == "filtered_read_mapping":
        filtered_read_mapping(exp)
    elif args.tool == "personalized_read_mapping":
        personalized_read_mapping(exp)
    elif args.tool == "create_multitarget_vcfs":
        create_multitarget_vcfs(exp, args.target_exp)
    elif args.tool == "VCFtoNP_parallel":
        VCFtoNP_parallel(exp)

if __name__ == '__main__':
    main()
