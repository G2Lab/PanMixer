using CategoricalArrays, Distances, Printf

"""
Function that implements the LD decay metric for ABC
"""

function get_closest_bp_to_cm(sorted_positions, bp_to_cm_map, bp::Int)::Float64
    # Find the closest position using binary search
    idx = searchsortedlast(sorted_positions, bp)

    if idx == 0
        # If the position is below the smallest key, extrapolate below
        closest_bp = sorted_positions[1]
        return bp_to_cm_map[closest_bp] - (sorted_positions[1] - bp) / 1_000_000
    elseif idx == length(sorted_positions)
        # If the position is above the largest key, extrapolate above
        closest_bp = sorted_positions[end]
        return bp_to_cm_map[closest_bp] + (bp - sorted_positions[end]) / 1_000_000
    else
        # If the position is between two keys, check for exact match or interpolate
        lower_bp = sorted_positions[idx]
        upper_bp = sorted_positions[idx + 1]

        if bp == lower_bp
            return bp_to_cm_map[lower_bp]  # Exact match
        elseif bp == upper_bp
            return bp_to_cm_map[upper_bp]  # Exact match
        elseif bp < upper_bp
            # Interpolate between lower_bp and upper_bp
            diff_bp = bp - lower_bp
            return bp_to_cm_map[lower_bp] + diff_bp / 1_000_000
        else
            # Extrapolate above
            diff_bp = bp - upper_bp
            return bp_to_cm_map[upper_bp] + diff_bp / 1_000_000
        end
    end
end

function filter_long_lines(input_path::String, output_path::String; max_length::Int = 1000)
    open(output_path, "w") do outfile
        open(input_path, "r") do infile
            for line in eachline(infile)
                if length(line) <= max_length
                    println(outfile, line)
                end
            end
        end
    end
end

function load_ref_LD_decay(cache_file::String)
    # Load the LD decay data from the reference file
    @info "Loading cached LD decay from $cache_file"
    df = DataFrame(CSV.File(cache_file))
    return [df.dist df.R2_mean]
end

function LD_decay(plink_file, plink, mapthin, out_prefix, bp_to_cm_map; cache_file=nothing)
    # Check for cache
    #if cache_file !== nothing && isfile(cache_file)
    #    @info "Loading cached LD decay from $cache_file"
    #    df = DataFrame(CSV.File(cache_file))
    #    return [df.dist df.R2_mean]
    #end

    # thin snp list to speed up calculations
    bim_file = @sprintf("%s.bim", plink_file)
    snp_thin = @sprintf("%s_snp_thin", out_prefix)
    run(`$mapthin -b 100 $bim_file $snp_thin`)

    snp_thin_filtered = @sprintf("%s_snp_thin_filtered", out_prefix)
    filter_long_lines(snp_thin, snp_thin_filtered; max_length=1000)

    # compute r2 values between all pairs
    run(`$plink --bfile $plink_file --extract $snp_thin_filtered --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out $snp_thin_filtered --allow-extra-chr`)
    
    # load LD data
    ld_df = DataFrame(CSV.File(@sprintf("%s.ld", snp_thin_filtered), delim=' ', ignorerepeated=true))

    # filter variants based on map
    total_variants = nrow(ld_df)  # Total number of variants before filtering

    @info "Filtering variants based on map"
    sorted_positions = sort(collect(keys(bp_to_cm_map)))
    @info "Variant positions sorted"
    ld_df.CM_A = [get_closest_bp_to_cm(sorted_positions, bp_to_cm_map, x) for x in ld_df.BP_A]
    ld_df.CM_B = [get_closest_bp_to_cm(sorted_positions, bp_to_cm_map, x) for x in ld_df.BP_B]

    # compute genetic distance
    ld_df.dist = ld_df.CM_B .- ld_df.CM_A

    # bin into intervals
    interval = 0.01
    ld_df.bin = cut(ld_df.dist, Vector(0:interval:maximum(ld_df.dist)), extend=missing)

    # group and average
    gdf = groupby(ld_df, :bin)
    ld_grp = combine(gdf, :R2 => mean)
    ld_grp.dist = Vector(0:interval:interval * (nrow(ld_grp)-1))

    CSV.write(cache_file, ld_grp)
    @info "Cached LD decay to $cache_file"

    return [ld_grp.dist[1:end] ld_grp.R2_mean[1:end]]
end


function run_ld_decay(synfile, realfile, plink, mapthin, out_prefix, bp_to_cm_map)
    real_cache = joinpath(dirname(realfile), "real_ld_decay.csv")
    syn_cache = joinpath(dirname(synfile), "syn_ld_decay.csv")

    ld_decay_real = load_ref_LD_decay(real_cache)
    ld_decay_syn = LD_decay(synfile, plink, mapthin, out_prefix * "_syn", bp_to_cm_map; cache_file=syn_cache)

    min_size = min(size(ld_decay_real, 1), size(ld_decay_syn, 1))
    min_size = min(min_size, 100)  # Ensure max 100 points
    ld_decay_real = ld_decay_real[1:min_size, :]
    ld_decay_syn = ld_decay_syn[1:min_size, :]

    df = DataFrame(real_x=ld_decay_real[:,1], real_y=ld_decay_real[:,2], syn_x=ld_decay_syn[:,1], syn_y=ld_decay_syn[:,2])
    outfile = joinpath(dirname(out_prefix), "results-ld-decay.csv")
    CSV.write(outfile, df)
    @info "LD decay data saved at $outfile"

    outfile = joinpath(dirname(out_prefix), "results-ld-decay.png")
    fig = Plots.plot(size=(400, 400))
    plot!(fig, ld_decay_real[:,1], ld_decay_real[:,2], label="real", xaxis="Genetic distance (cm)", yaxis="LD estimate (r2)", title="LD decay")
    plot!(fig, ld_decay_syn[:,1], ld_decay_syn[:,2], label="synthetic")
    savefig(fig, outfile)
    @info "LD decay plot saved at $outfile"

    ld_distance = evaluate(Euclidean(), ld_decay_real[:,1], ld_decay_syn[:,1])
    @info @sprintf("LD decay distance between real and synthetic data is %f", ld_distance)
end
