"""Executes the pipeline for evaluating synthetic data quality
"""

include("metrics/eval_ld_decay.jl")
include("metrics/eval_maf.jl")

using Printf

function run_ld_decay_evaluation(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map)
    run_ld_decay(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map)
end


function run_maf_evaluation(real_maf_file, synt_maf_file)
    run_maf(real_maf_file, synt_maf_file)
end


function create_bp_cm_ref(genetic_distfile)
    # Read headerless tab-delimited file and assign column names
    col_names = [:rsid, :position, :genetic_position]
    bp_to_cm_df = CSV.read(genetic_distfile, DataFrame; delim=' ', header=false)
    # print column names
    rename!(bp_to_cm_df, col_names)

    # Create mapping from base pair position to genetic position
    bp_to_cm_map = Dict(zip(bp_to_cm_df.position, bp_to_cm_df.genetic_position))
    return bp_to_cm_map
end



"""Computations for MAF using PLINK
"""
function run_maf_tools(plink, reffile_prefix, synfile_prefix, ref_outdir, syn_outdir)
    @info "Running external tools for MAF"
    reffile_out = @sprintf("%s.ref.maf", ref_outdir)
    synfile_out = @sprintf("%s.syn.maf", syn_outdir)

    run(`$plink --bfile $synfile_prefix --freq --out $synfile_out -allow-extra-chr`)
    real_maffile = @sprintf("%s.frq", reffile_out)
    syn_maffile = @sprintf("%s.frq", synfile_out)
    return real_maffile, syn_maffile
end

"""Compute bfiles
"""
function run_bfiles(plink, reffile_prefix, synfile_prefix, ref_outdir, syn_outdir)
    @info "Running external tools for bfiles"
    reffile_out = @sprintf("%s.ref", ref_outdir)
    synfile_out = @sprintf("%s.syn", syn_outdir)

    run(`$plink --vcf $synfile_prefix --make-bed --out $synfile_out --allow-extra-chr --vcf-half-call missing`)
    return reffile_out, synfile_out
end


"""Entry point to running the evaluation pipeline for genotype data

Note that the evaluation pipeline assumes that the synthetic data you want 
to evaluate has already been geneerated, using the setup specified in the
configuration file. It is therefore recommended to run the pipeline with the 
--genotype and --evaluation flags together, so that the program generates 
the data and then immediately evaluates it using the correct settings.
"""
function run_evaluation(experiment_number, chromosome, row_id)
    plink_path = "/gpfs/commons/home/jblindenbach/tools/plink"
    mapthin_path = "/gpfs/commons/home/jblindenbach/tools/mapthin-v1.11-linux-x86_64/mapthin"

    reffile_prefix = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr$chromosome/pangenome.vcf.gz"
    synfile_prefix = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_$experiment_number/data/chr$chromosome/$row_id/new_haplotypes.vcf.gz"

    # Check if the files exist
    if !isfile(synfile_prefix)
        synfile_prefix = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_$experiment_number/data/chr$chromosome/$row_id/new_haplotypes.vcf"
        if !isfile(synfile_prefix)
            @error "Synthetic file not found: $synfile_prefix"
            return
        end
    end

    outdir = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_$experiment_number/data/chr$chromosome/$row_id/"
    outdir_ref = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr$chromosome/pangenome"
    outdir_syn = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_$experiment_number/data/chr$chromosome/$row_id/new_haplotypes"

    reffile_bfile, synfile_bfile = run_bfiles(plink_path, reffile_prefix, synfile_prefix, outdir_ref, outdir_syn)

    reffile_maf, synfile_maf = run_maf_tools(plink_path, reffile_bfile, synfile_bfile, outdir_ref, outdir_syn)
    run_maf_evaluation(reffile_maf, synfile_maf)

    bp_to_cm_map = create_bp_cm_ref("/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/1000-genomes-genetic-maps/interpolated_from_hapmap/chr$chromosome.interpolated_genetic_map")

    run_ld_decay_evaluation(synfile_bfile, reffile_bfile, plink_path, mapthin_path, outdir, bp_to_cm_map)

    clean_up(outdir)
end

function run_setup_for_reference(chromosome)
    plink_path = "/gpfs/commons/home/jblindenbach/tools/plink"
    mapthin_path = "/gpfs/commons/home/jblindenbach/tools/mapthin-v1.11-linux-x86_64/mapthin"

    reffile_prefix = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr$chromosome/pangenome.vcf.gz"
    outdir_ref = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr$chromosome/pangenome"

    reffile_out = @sprintf("%s.ref", outdir_ref)
    reffile_out_freq = @sprintf("%s.ref.maf", outdir_ref)

    run(`$plink_path --vcf $reffile_prefix --make-bed --out $reffile_out --allow-extra-chr --vcf-half-call missing`)
    run(`$plink_path --bfile $reffile_out --freq --out $reffile_out_freq --allow-extra-chr`)

    bp_to_cm_map = create_bp_cm_ref("/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/1000-genomes-genetic-maps/interpolated_from_hapmap/chr$chromosome.interpolated_genetic_map")
    
    real_cache = joinpath(dirname(reffile_out), "real_ld_decay.csv")

    ld_decay_real = LD_decay(reffile_out, plink_path, mapthin_path, outdir_ref * "_real", bp_to_cm_map; cache_file=real_cache)
end



function clean_up(dir::AbstractString)
    @info "Cleaning up temporary files in directory: $dir"

    # List of extensions to delete
    extensions = [".bed", ".bim", ".fam", ".log", ".nosex", ".raw", ".frq", ".ld"]

    # Scan the directory for files with matching extensions
    for file in readdir(dir; join=true)
        for ext in extensions
            if endswith(file, ext)
                rm(file; force=true)
                @info "Deleted: $file"
                break
            end
        end
    end

    return nothing
end

function main()
    # Get command-line arguments
    args = ARGS
    if length(args) < 2
        println("Usage: julia script.jl <experiment_number> <chromosome> <row_id>")
        return
    end

    if args[1] == "ref"
        chromosome = args[2]
        run_setup_for_reference(chromosome)
        return
    end

    experiment_number = parse(Int, args[1])
    chromosome =args[2]
    row_id = args[3]


    run_evaluation(experiment_number, chromosome, row_id)
end

# Call main if this script is run directly
main()