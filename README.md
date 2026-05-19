# PanMixer Toolkit

A command-line toolkit for running privacy–utility experiments on pangenome graphs.  
This pipeline supports experiment initialization, optimization, stacking, VCF merging, and a wide set of downstream analyses.

This toolkit is designed to run in a cluster environment and currently supports execution via `slurm` using `sbatch`.

> This repository accompanies the PanMixer paper. The `main` branch reflects the publication-ready state of the toolkit.

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
conda create -y -n panmixer python=3.11
conda activate panmixer
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

We suggest upon cloning this repository to copy the config, `cp config.yaml config.local.yaml` and changing the configurations in `config.local.yaml`

> Note: Large data assets (pangenomes, read fastqs, 1000 Genomes datasets) are **not bundled** in this repository. You must download them yourselves. We provide scripts to help you :)
---

## Dataset installation and preprocessing

`cd` into the starting-data directory `./starting_data/`. Run each command in `get_data.sh` one at a time, waiting for the previous command to finish.

> WARNING: External download links may break over time. If a script fails, check the source URL inside it before retrying.

Commands:
- `./scripts/get_pangenomes.sh`: Downloads the PGGB draft human pangenome
- `./scripts/get_pangenie_alignments.sh`: Downloads the Pangenie callset
- `./scripts/get_1000g_phased.sh`: Downloads the 1000 Genomes Project phased panel
- `bcftools index pangenome.vcf.gz`: Indexes the downloaded pangenome
- `./scripts/remove_X.sh`: Removes any X chromosome data
- `./scripts/remove_chm13.sh`: Removes the chm13 backbone
- `sbatch scripts/split_data.sbatch`: Splits the data by chromosome
- `sbatch scripts/get_num_alleles.sbatch`: Identifies the unique alleles and variants
- `sbatch scripts/get_blocks.sbatch`: Computes the LD blocks
- `sbatch scripts/convert_2_npy.sbatch`: Converts the VCF files into Numpy files for easier IO
- `sbatch scripts/get_mappings.sbatch`: Computes variant mappings
- `sbatch scripts/segment_blocks.sbatch`: Refines segmented blocks produced by plink
- `sbatch scripts/get_pmi_utility.sbatch`: Computes the PMI and utility loss for each obfuscation move

These commands take about 1–2 days to run end-to-end. Please be patient.

To make the sampling step reproducible during preprocessing, pass a base seed to the SLURM job:

```bash
sbatch --export=SEED=123 scripts/get_pmi_utility.sbatch
```

Each array task derives its own seed from this base seed, so different chromosome/subject tasks do not reuse the same random stream.

## Quick Start

```bash
# 1) Start a new experiment
python3 main.py experiment_starter
# 2) Run optimizer
python3 main.py --exp 0 optimize --fixed_param utility

# Optional: make random optimizer baselines reproducible
python3 main.py --exp 0 --seed 123 optimize --fixed_param random

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
- `--seed INT` – base random seed for commands that sample. Currently used by the random optimizer baseline; per-SLURM-task seeds are derived from this base seed.

### Reproducibility

The HMM and allele-frequency resampling paths support explicit seeds. Without a seed, NumPy uses its default process-level randomness.

For the standalone obfuscation tool:

```bash
python3 tools/obfuscate.py starting_data HG00438 0.1 21 obfuscation_output --seed 123
```

For the obfuscation SLURM wrapper:

```bash
sbatch --export=DATA_DIR=/path/to/starting_data,SUBJECT=HG00438,CAPACITY=0.1,OUTPUT_DIR=/path/to/output,SEED=123 tools/obfuscate.sbatch
```

For starting-data PMI/utility precomputation:

```bash
sbatch --export=SEED=123 starting_data/scripts/get_pmi_utility.sbatch
```

### Subcommands

#### **experiment_starter**
Initialize an experiment.

```bash
python3 main.py experiment_starter \
  --capacity_file path/to/capacities.txt \
  --subjects_file path/to/subjects.txt \
  [--baseline_unedited] [--baseline_empty] [--baseline_unique]
```

Parameters:
- `--capacity_file` supplies a new line delimited file with all the target capacities (target privacy risk and utility loss)
- `--subjects_file` supplies a list of target individuals to obfuscate
- `--baseline_unedited` runs an experiment with the original pangenome graphs (no edits) as a baseline
- `--baseline_empty` runs an experiment with the subject removed
- `--baseline_unique` runs a baseline experiment where unique-to-subject variants are removed

#### **optimize**
Run linear optimizer.

```bash
python3 main.py optimize --fixed_param utility
```

Parameters:
- `--fixed_param` can be either `privacy` or `utility` and optimize for a fixed capacity limit as set by the capacity file of that parameter
- `--baseline_unique` optimize for the baseline with removed unique variants (use with experiments started with `--baseline_unique`)

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

#### **af_loss / ld_loss**
Compute allele frequency loss and linkage disequilibrium loss.

```bash
python3 main.py af_loss
python3 main.py ld_loss
```

Parameters (both commands):
- `--dont_replace` compute metrics without replacing existing results

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



#### **verify**
Verify genotype consistency across all samples for an experiment (runs per-chromosome in parallel via Slurm).

```bash
python3 main.py --exp 0 verify
```

#### **VCFtoNP_parallel**
Convert per-subject VCF files to NumPy arrays in parallel across chromosomes via Slurm. Useful for preprocessing new VCF datasets into the internal numpy format used by the pipeline.

```bash
python3 main.py --exp 0 VCFtoNP_parallel
```

#### **create_multitarget_vcfs**
Create multi-target VCFs by merging obfuscation results from a source experiment into the format expected by a target experiment. Useful for composing experiments that cover multiple subjects.

```bash
python3 main.py --exp 0 create_multitarget_vcfs --target_exp 1
```

Parameters:
- `--target_exp` the experiment number to use as the target layout

#### **gap_threshold**
Analyze gap score privacy thresholds across a range of pangenome sizes. Useful for understanding how privacy guarantees scale with pangenome size.

```bash
python3 main.py --exp 0 gap_threshold --pangenome_sizes 50,100,250,500,1000,2000
```

Parameters:
- `--pangenome_sizes` comma-separated list of pangenome sizes to test (default: `50,100,250,500,1000,1500,2000,2500,3000`)

#### **gather_results**
Aggregate outputs across experiments.

```bash
python3 main.py [--overwrite] gather_results \
  [--optimizer] [--reindex] [--gap_score] [--stacker] \
  [--af_loss] [--ld_loss] [--pangenie_stats] \
  [--accuracy_stats] [--giraffe]
```

Run without flags to gather all results, or pass one or more flags to gather only those result types.

Parameters:
- `--overwrite` overwrites the results if the flag is present (reset)
- `--optimizer` gather optimizer results
- `--reindex` re-index the experiment
- `--gap_score` gather gap score results
- `--stacker` gather stacker results
- `--af_loss` gather allele frequency loss results
- `--ld_loss` gather LD loss results
- `--pangenie_stats` gather Pangenie stats
- `--accuracy_stats` gather Beagle accuracy stats
- `--giraffe` gather Giraffe alignment results
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
   ├─ gap_threshold
   ├─ af_loss / ld_loss
   ├─ beagle / accuracy_stats
   ├─ vg_prep / quick_align
   ├─ verify
   └─ gather_results

[ Utilities ]
   ├─ VCFtoNP_parallel  (preprocess new VCF data)
   └─ create_multitarget_vcfs  (compose multi-subject experiments)
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
