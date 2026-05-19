"""
Combine per-chromosome obfuscated VCFs into a single VCF and aggregate stats.

Usage:
    python tools/combine_obfuscated.py <output_dir>

Expects output_dir/chr{1..22}/obfuscated.vcf.gz from obfuscate.py.
Produces output_dir/obfuscated_all.vcf.gz and output_dir/stats_all.json.
"""

import json
import os
import sys
import subprocess


def main():
    output_dir = sys.argv[1]

    # Collect per-chromosome VCFs
    vcfs = []
    stats_all = []
    for c in range(1, 23):
        vcf = f"{output_dir}/chr{c}/obfuscated.vcf.gz"
        stat = f"{output_dir}/chr{c}/stats.json"
        if not os.path.exists(vcf):
            print(f"WARNING: missing {vcf}")
            continue
        vcfs.append(vcf)
        if os.path.exists(stat):
            stats_all.append(json.load(open(stat)))

    if not vcfs:
        print("No chromosome VCFs found.")
        sys.exit(1)

    # Concatenate VCFs
    out_vcf = f"{output_dir}/obfuscated_all.vcf.gz"
    cmd = ["bcftools", "concat", "-Oz", "-o", out_vcf] + vcfs
    print(f"Concatenating {len(vcfs)} VCFs...", flush=True)
    subprocess.run(cmd, check=True)
    subprocess.run(["bcftools", "index", out_vcf], check=True)

    # Aggregate stats
    if stats_all:
        agg = {
            "subject_name": stats_all[0]["subject_name"],
            "capacity": stats_all[0]["capacity"],
            "chromosomes": len(stats_all),
            "total_utility_loss": sum(s["utility_loss"] for s in stats_all),
            "total_max_utility_loss": sum(s["max_utility_loss"] for s in stats_all),
            "total_pmi_gain": sum(s["pmi_gain"] for s in stats_all),
            "total_max_pmi_gain": sum(s["max_pmi_gain"] for s in stats_all),
            "total_moves": sum(s["moves"] for s in stats_all),
            "total_changed_alleles": sum(s["changed_alleles"] for s in stats_all),
            "total_alleles": sum(s["total_alleles"] for s in stats_all),
        }
        agg["utility_loss_frac"] = agg["total_utility_loss"] / agg["total_max_utility_loss"]
        agg["pmi_gain_frac"] = (
            agg["total_pmi_gain"] / agg["total_max_pmi_gain"]
            if agg["total_max_pmi_gain"] > 0 else 0
        )
        agg["changed_frac"] = agg["total_changed_alleles"] / agg["total_alleles"]

        # Aggregate timings and RSS across chromosomes
        timing_keys = ["load_data_s", "graph_prep_s", "optimizer_s", "stacker_s", "vcf_write_s", "total_s"]
        rss_keys = ["load_data_rss_mb", "graph_prep_rss_mb", "optimizer_rss_mb", "stacker_rss_mb", "vcf_write_rss_mb", "peak_rss_mb"]

        timings_agg = {}
        for key in timing_keys:
            vals = [s["timings"][key] for s in stats_all if "timings" in s and key in s["timings"]]
            if vals:
                timings_agg[f"sum_{key}"] = round(sum(vals), 2)
                timings_agg[f"max_{key}"] = round(max(vals), 2)
                timings_agg[f"avg_{key}"] = round(sum(vals) / len(vals), 2)

        for key in rss_keys:
            vals = [s["timings"][key] for s in stats_all if "timings" in s and key in s["timings"]]
            if vals:
                timings_agg[f"max_{key}"] = round(max(vals), 1)

        agg["timings"] = timings_agg
        agg["per_chromosome"] = [
            {
                "chromosome": s["chromosome"],
                "timings": s.get("timings", {}),
            }
            for s in stats_all
        ]

        json.dump(agg, open(f"{output_dir}/stats_all.json", "w"), indent=2)
        print(f"\nAggregated stats:")
        print(f"  Utility loss: {agg['utility_loss_frac']*100:.1f}%")
        print(f"  PMI gain: {agg['pmi_gain_frac']*100:.1f}%")
        print(f"  Changed alleles: {agg['changed_frac']*100:.2f}%")
        print(f"\nTimings (sum across {len(stats_all)} chromosomes):")
        for key in timing_keys:
            sk = f"sum_{key}"
            mk = f"max_{key}"
            if sk in timings_agg:
                print(f"  {key:<18} sum={timings_agg[sk]:>8.1f}s  max={timings_agg[mk]:>8.1f}s")
        print(f"\nPeak RSS (max across chromosomes):")
        for key in rss_keys:
            mk = f"max_{key}"
            if mk in timings_agg:
                print(f"  {key:<24} {timings_agg[mk]:>8.1f} MB")

    print(f"\nOutput: {out_vcf}")


if __name__ == "__main__":
    main()
