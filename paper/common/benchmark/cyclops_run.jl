# Driver for CYCLOPS 2.0 (Hammarlund & Anafi). Runs CYCLOPS.Fit and CYCLOPS.Align on a genes x samples CSV.
# Usage: julia --project=<cyclops_env> cyclops_run.jl expr=<csv> seeds=<txt> out=<dir> [hyper=<toml>] [times=<csv>]
#        [seed=0] [workers=8] [cyclops=<path to CYCLOPS.jl>]
# expr:  first column gene symbols, header = sample IDs (letters/digits/underscore only), non-log expression.
# seeds: one seed gene per line.   times: CSV with columns sample,time_hours (known collection times, any subset).
# hyper: TOML overriding the CYCLOPS 2.0 template hyperparameters (keys as in the template, without the leading colon).
# Outputs in <out>: phases.csv (sample, phase_raw, phase, aligned_by), fit_output.csv, eigengene_correlations.csv,
#                   run_info.json, cyclops_raw/ (CYCLOPS.Align's own Plots/Fits/Models/Parameters).
using Distributed, Random, LinearAlgebra, TOML, Dates
args = Dict(split(a, "="; limit = 2) |> x -> (x[1] => x[2]) for a in ARGS)
nw = parse(Int, get(args, "workers", "8"))
cyclops_path = get(args, "cyclops", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tools/CYCLOPS-2.0/CYCLOPS.jl")
ENV["MPLBACKEND"] = "Agg"   # non-interactive backend for the (unavoidable) Align plots
addprocs(nw; exeflags = ["--project=$(Base.active_project())"])
@everywhere ENV["MPLBACKEND"] = "Agg"
@everywhere using LinearAlgebra
@everywhere LinearAlgebra.BLAS.set_num_threads(1)
LinearAlgebra.BLAS.set_num_threads(nw)
@everywhere include($cyclops_path)
using DataFrames, CSV, Statistics

seed = parse(Int, get(args, "seed", "0"))
out = args["out"]; mkpath(out)
expr = CSV.read(args["expr"], DataFrame; stringtype = String)
seed_genes = String.(filter(!isempty, strip.(readlines(args["seeds"]))))

# Training hyperparameters: defaults of CYCLOPS_2_0_Template.jl.
tp = Dict{Symbol,Any}(:regex_cont => r".*_C", :regex_disc => r".*_D", :blunt_percent => 0.975,
    :seed_min_CV => 0.14, :seed_max_CV => 0.9, :seed_mth_Gene => 10000,
    :norm_gene_level => true, :norm_disc => false, :norm_disc_cov => 1,
    :eigen_reg => true, :eigen_reg_disc_cov => 1, :eigen_reg_exclude => false, :eigen_reg_r_squared_cutoff => 0.6,
    :eigen_reg_remove_correct => false, :eigen_first_var => false, :eigen_first_var_cutoff => 0.85,
    :eigen_total_var => 0.85, :eigen_contr_var => 0.05, :eigen_var_override => true, :eigen_max => 5,
    :out_covariates => true, :out_use_disc_cov => true, :out_all_disc_cov => true, :out_disc_cov => 1,
    :out_use_cont_cov => false, :out_all_cont_cov => true, :out_use_norm_cont_cov => false,
    :out_all_norm_cont_cov => true, :out_cont_cov => 1, :out_norm_cont_cov => 1,
    :init_scale_change => true, :init_scale_1 => false,
    :train_n_models => 80, :train_μA => 0.001, :train_β => (0.9, 0.999), :train_min_steps => 1500,
    :train_max_steps => 2050, :train_μA_scale_lim => 1000, :train_circular => false,
    :train_collection_times => true, :train_collection_time_balance => 1.0,
    :cosine_shift_iterations => 192, :cosine_covariate_offset => true,
    :align_p_cutoff => 0.05, :align_base => "radians", :align_disc => false, :align_disc_cov => 1,
    :align_other_covariates => false, :align_batch_only => false,
    :X_Val_k => 10, :X_Val_omit_size => 0.1,
    :plot_use_o_cov => true, :plot_correct_batches => true, :plot_disc => false, :plot_disc_cov => 1,
    :plot_separate => false, :plot_color => ["b", "orange", "g", "r", "m", "y", "k"],
    :plot_only_color => true, :plot_p_cutoff => 0.05)
for (k, v) in (haskey(args, "hyper") ? TOML.parsefile(args["hyper"]) : Dict())
    tp[Symbol(k)] = k == "train_β" ? Tuple(Float64.(v)) : v
end
# CYCLOPS.GeneLevelCutoffs indexes the (seed_mth_Gene+1)-th most expressed gene: needs more than that many genes.
ngenes = nrow(expr)
tp[:seed_mth_Gene] = min(tp[:seed_mth_Gene], ngenes - 1)

# Time-stamped samples used during training (semi-supervised; CYCLOPS wants radians).
if haskey(args, "train_times")
    tr = CSV.read(args["train_times"], DataFrame; stringtype = String)
    tp[:train_sample_id] = String.(tr[:, 1]); tp[:train_sample_phase] = mod.(Float64.(tr[:, 2]) .* (2π / 24), 2π)
end

# Known collection times (hours) -> post-training alignment (CYCLOPS.Align, :align_samples / :align_phases).
if haskey(args, "times")
    tt = CSV.read(args["times"], DataFrame; stringtype = String)
    tp[:align_samples] = String.(tt[:, 1]); tp[:align_phases] = Float64.(tt[:, 2]); tp[:align_base] = "hours"
end

# CYCLOPS.Fit re-seeds Random with 1234 internally; the seed here instead fixes the sample (column) order,
# which changes the sequence of per-sample Adam updates and therefore the trained model.
Random.seed!(seed)
ids = names(expr)[2:end]
perm = randperm(length(ids))
expr = expr[:, vcat(1, 1 .+ perm)]

t0 = time()
eigendata, fit, corr, model, ops = CYCLOPS.Fit(expr, seed_genes, tp)
CYCLOPS.Align(expr, fit, corr, model, ops, joinpath(out, "cyclops_raw"))
runtime = time() - t0

fit[!, :Phase] = Float64.(fit.Phase)
aligned_by = hasproperty(fit, :Phases_SA) ? "samples" : hasproperty(fit, :Phases_MA) ? "clock_genes" : "none"
phase = aligned_by == "samples" ? fit.Phases_SA : aligned_by == "clock_genes" ? fit.Phases_MA : fit.Phase
res = DataFrame(sample = fit.ID, phase_raw = fit.Phase, phase = mod.(Float64.(phase), 2π), aligned_by = aligned_by)
CSV.write(joinpath(out, "phases.csv"), res)
CSV.write(joinpath(out, "fit_output.csv"), fit)
CSV.write(joinpath(out, "eigengene_correlations.csv"), corr)
open(joinpath(out, "run_info.json"), "w") do io
    println(io, "{\"runtime_s\": $runtime, \"n_seed_genes_used\": $(length(ops[:o_seed_genes])), ",
        "\"n_eigengenes\": $(ops[:o_svd_n_dims]), \"cumvar\": $(round.(ops[:o_svd_cumvar][1:min(10,end)], digits=4)), ",
        "\"aligned_by\": \"$aligned_by\", \"seed\": $seed, \"julia\": \"$VERSION\"}")
end
