# PanMixer Toolkit

A command-line toolkit for running privacy–utility experiments on pangenome graphs.  
This pipeline supports experiment initialization, optimization, stacking, VCF merging, and a wide set of downstream analyses.

This toolkit is designed to run in a cluster environment and currently supports execution via `slurm` using `sbatch`.

---

## Overview

The toolkit organizes workflows around **experiments** identified by an `--exp` number.  
If not specified, the latest experiment is automatically selected.

**Typical flow:**
1. `experiment_starter` 
2. `optimize` 
3. `stacker`
4. `convert_2_vcf`
5. `combine_vcfs`
6. **Downstream analyses** (`gap_score`, `af_loss`, `ld_loss`, `beagle`, etc.)
7. `gather_results`

---

## Installation

```bash
# Create and activate an environment
conda create -y -n panmixer4 python=3.11
conda activate panmixer4
# or use venv:
# python -m venv .venv && source .venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt

# Ensure external tools (bcftools, tabix, plink) are available in PATH
```

## Data & Paths

Defined in `config.yaml`, these constants need to be set before running any scripts
- `python_env` – path to python environment (ex: ./envs/panmixer-env)
- `base_path` – base path to root of this directory (ex: ./PanMixer)

> Note: Large data assets (pangenomes, read fastqs, 1000 Genomes datasets) are **not bundled** in this repository. You must download them yourselves. We provide scripts to help you :)
---

## Dataset installation and preprocessing



## Quick Start

```bash
# 1) Start a new experiment
python3 main.py experiment_starter
# 2) Run optimizer
python3 main.py --exp 0 optimize --fixed_param utility

# 3) Stack edits using a strategy
python3 main.py --exp 0 stacker --strategy to_best

# 4) Merge per-chromosome VCFs into one
python3 main.py --exp 0 convert_2_vcf
python3 main.py --exp 0 combine_vcfs

# 5) Run downstream analyses
python3 main.py --exp 0 gap_score_all
python3 main.py --exp 0 af_loss
python3 main.py --exp 0 ld_loss
python3 main.py --exp 0 vg_prep
python3 main.py --exp 0 quick_align
```

## Command Reference

### Global Flags
- `--exp INT` – experiment number (defaults to latest if not set).
- `--overwrite` – allow overwriting existing results.

### Subcommands

#### **experiment_starter**
Initialize an experiment.

```bash
python3 main.py experiment_starter   --capacity_file path/to/capacities.txt   --subjects_file path/to/subjects.txt   [--baseline_unedited] [--baseline_empty]
```

Parameters:
- `--capacity_file` supplies a new line delimited file with all the target capacities (target privacy risk and utility loss)
- `--subjects_file` supplies a list of target individuals to obfuscate
- `--baseline_unedited` runs an experiment with the original pangenome graphs (no edits) as a baseline
- `--baseline_empty` runs an experiment with the subject removed


#### **optimize**
Run linear optimizer.

```bash
python3 main.py optimize --fixed_param utility
```

Parameters:
- `--fixed_param` can be either `privacy` or `utility` and optimize for a fixed capacity limit as set by the capacity file of that parameter

#### **stacker**
Apply stacking strategy which applies the obfuscated moves taken by the optimizer step.

```bash
python3 main.py stacker --strategy to_best
```

Parameters:
- `--strategy` can either be `to_best` which is the strategy described in the paper, `to_empty` removes the individual and is run with experiments setup using the `--baseline_empty` flag, or `to_unedited` which ignores all obfuscation moves and is used with `--baseline_unedited`

#### **convert_2_vcf**
Converts the output of stacker to a variant centric representation of the graph
```bash
python3 main.py convert_2_vcf
```

#### **combine_vcfs**
Merge per-chromosome VCFs.

```bash
python3 main.py combine_vcfs
```

#### **gap_score / gap_score_all**
Compute haploid/diploid gap scores.

```bash
python3 main.py gap_score
python3 main.py gap_score_all
```

Run `gap_score` for one attack
Run `gap_score_all` for all attacks

#### **af_loss / ld_loss / maf_ld**
Compute allele frequency loss, linkage disequilibrium loss.

```bash
python3 main.py af_loss
python3 main.py ld_loss
```

#### **Beagle**
Computes beagle reconstruction of genotypes
```bash
python3 main.py beagle
```

#### **accuracy_stats**
Get accuracy statistics of beagle refinements.

```bash
python3 main.py accuracy_stats
```

#### **vg_prep / quick_align**
Prepare `vg giraffe` indices and run read mapping for all reads in `/read_fastqs/`.

```bash
python3 main.py vg_prep
python3 main.py quick_align
```

#### **gather_results**
Aggregate outputs across experiments.

```bash
python3 main.py [--overwrite] gather_results [--optimizer] [--reindex] [--gap_score] [--stacker]   [--af_loss] [--ld_loss] [--pangenie_stats]   [--accuracy_stats] [--maf_ld] [--giraffe]   
```

Gets results from all of the above commands, can either be run with the flag of the command or without to get all results

Parameters:
- `--overwrite` overwrites the results if the flag is present (reset)
---

## Workflow Diagram

```
[ experiment_starter ]
        ↓
     [ optimize ]
        ↓
      [ stacker ]
        ↓
    [ convert_2_vcf ]
        ↓
   [ combine_vcfs ]
        ↓
[ downstream analyses ]
   ├─ gap_score / gap_score_all
   ├─ af_loss / ld_loss
   ├─ beagle / beagle_stats / accuracy_stats
   ├─ vg_prep / quick_align
   └─ gather_results
```

---

## Outputs

- **Experiments**: stored in directories per experiment number in `/experiments/`.
- **VCFs**: per-chromosome and merged VCFs.
- **Metrics**: AF/LD loss, gap scores.
- **Alignments**: vg-prepped indexes, alignment results.
- **Aggregated results**: via `gather_results`.
---

## Troubleshooting

- **No experiment found**: Run `experiment_starter` first.
- **Missing input files**: Verify `STARTING_DATA_PATH`, `--subjects_file`, and `--capacity_file`.
- **VCF issues**: Ensure inputs are bgzipped (`.vcf.gz`) and indexed (`.tbi`).
- **External tool errors**: Ensure `bcftools`, `tabix`, `plink` are installed and on PATH.
---